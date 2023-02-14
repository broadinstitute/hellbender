import argparse
import os
import re

from azure.identity import DefaultAzureCredential
from azure.mgmt.resource import ResourceManagementClient
from azure.mgmt.storage import StorageManagementClient
from azure.mgmt.subscription import SubscriptionClient
from azure.storage.blob import BlobClient, BlobServiceClient

from pathlib import Path
from uuid import uuid4


def exactly_one_or_die(result, kind, filter=None, describer=None):
    found = list(result)

    if filter:
        found = [f for f in found if filter(f)]

    if len(found) == 0:
        raise ValueError(f"Could not find a {kind}!")

    if len(found) > 1:
        message = f"Found multiple {kind}!"
        if describer:
            message = message + " : " + ", ".join([describer(f) for f in found])
        raise ValueError(message)

    return found[0]


def get_subscription(credentials):
    subscription_client = SubscriptionClient(credentials)
    return exactly_one_or_die(subscription_client.subscriptions.list(),
                              "subscription",
                              describer=lambda s: s.id)


def get_resource_group(credentials, subscription):
    resource_client = ResourceManagementClient(credentials, subscription.subscription_id)

    if args.resource_group:
        resource_group_filter = lambda g: g.name == args.resource_group
        resource_group_descriptor = f"resource group '{args.resource_group}'"
    else:
        pattern = f"{os.environ['USER']}-[a-f0-9]+$"
        resource_group_descriptor = f"resource group matching pattern '{pattern}'"
        resource_group_filter = lambda g: re.match(pattern, g.name)

    return exactly_one_or_die(resource_client.resource_groups.list(),
                              resource_group_descriptor,
                              filter=resource_group_filter)


def get_storage_account(credentials, subscription, resource_group):
    storage_client = StorageManagementClient(credentials, subscription.subscription_id)
    # `az` returns storage account JSONs with a `resourceGroup` attribute, but the objects returned by the Python API do
    # not have this attribute. However the `id`s of these Python objects do contain the resource group in a predictable
    # pattern, so look for that instead.
    id_prefix = f"/subscriptions/{subscription.subscription_id}/resourceGroups/{resource_group.name}"
    return exactly_one_or_die(storage_client.storage_accounts.list(), "storage account",
                              filter=lambda a: a.id.startswith(id_prefix),
                              describer=lambda a: a.name)


def write_input_file(storage_account, path, data):
    blob_url = f"{storage_account.primary_endpoints.blob}/inputs/{path}"

    blob_client = BlobClient.from_blob_url(
        blob_url=blob_url,
        credential=credentials
    )

    blob_client.upload_blob(data)


def get_blob_service_client():
    return BlobServiceClient.from_connection_string(os.getenv('AZURE_CONNECTION_STRING'))


def generate_trigger_json(workflow_storage_path, inputs_storage_path):
    # If defined, enclose in double quotes for interpolation into the f string below.
    inputs_storage_path = f'"{inputs_storage_path}"' if inputs_storage_path else "null"

    return f"""
{{
  "WorkflowUrl": "{workflow_storage_path}",
  "WorkflowInputsUrl": {inputs_storage_path},
  "WorkflowInputsUrls": null,
  "WorkflowOptionsUrl": null,
  "WorkflowDependenciesUrl": null
}}
    """.strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Submit workflow to Cromwell on Azure')
    parser.add_argument('--workflow', type=str, help='Workflow WDL source', required=True)
    parser.add_argument('--inputs', type=str, help='Workflow inputs', required=False)
    parser.add_argument('--resource-group', type=str, help='Azure Resource Group name', required=False)
    args = parser.parse_args()

    if not os.getenv('AZURE_CONNECTION_STRING'):
        raise ValueError("Must define 'AZURE_CONNECTION_STRING' as a SAS token, see https://learn.microsoft.com/en-us/azure/storage/common/storage-configure-connection-string#store-a-connection-string")

    # https://github.com/Azure/azure-sdk-for-python/issues/22822#issuecomment-1024668507
    credentials = DefaultAzureCredential(exclude_shared_token_cache_credential=True)

    subscription = get_subscription(credentials)
    resource_group = get_resource_group(credentials, subscription)
    storage_account = get_storage_account(credentials, subscription, resource_group)
    blob_service_client = get_blob_service_client()

    inputs_client = blob_service_client.get_container_client('inputs')
    # 'name' is the filename without leading directory components.
    # e.g. name for /path/to/hello.wdl is hello.wdl
    # 'stem' is the filename without leading directory components and without an extension.
    # e.g. stem for /path/to/hello.wdl is hello
    workflow_path = Path(args.workflow)

    with open(args.workflow, "rb") as workflow:
        blob_client = inputs_client.get_blob_client(f"{workflow_path.stem}/{workflow_path.name}")
        blob_client.upload_blob(workflow)
    workflow_storage_path = f'/{storage_account.name}/inputs/{workflow_path.stem}/{workflow_path.name}'

    inputs_storage_path = None
    if args.inputs:
        inputs_path = Path(args.inputs)
        with open(args.inputs, "rb") as inputs:
            blob_client = inputs_client.get_blob_client(f"{workflow_path.stem}/{inputs_path.name}")
            blob_client.upload_blob(inputs)
        inputs_storage_path = f'/{storage_account.name}/inputs/{workflow_path.stem}/{workflow_path.name}'

    # Create the trigger JSON and stage into /storage account/workflows/new with a name of
    # {workflow_path.stem}-<random UUID>.json
    trigger_json = generate_trigger_json(workflow_storage_path, inputs_storage_path)
    workflows_client = blob_service_client.get_container_client('workflows')

    trigger_file_path = f'new/{workflow_path.stem}-{uuid4()}.json'
    blob_client = workflows_client.get_blob_client(trigger_file_path)
    blob_client.upload_blob(bytes(trigger_json, 'utf8'))

    print(f"Trigger JSON staged to /{storage_account.name}/workflows/{trigger_file_path}.")

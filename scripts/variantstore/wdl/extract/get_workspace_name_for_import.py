import json
import requests
import argparse

from terra_notebook_utils import table
from terra_notebook_utils import gs
# from terra_notebook_utils import workspace
import re

# make a default, but allow a user to overwrite it


def get_workspace_name(workspace_id, workspace_name_output, workspace_namespace_output):
    with open(workspace_name_output, "w") as name_output, open(workspace_namespace_output, "w") as namespace_output:
        token = gs.get_access_token()
        # grab the workspace information from rawls
        rawls = 'https://rawls.dsde-prod.broadinstitute.org/api/workspaces/id/{}?fields=workspace.namespace,workspace.googleProject'.format(workspace_id)
        head = {'Authorization': 'Bearer {}'.format(token)}
        response = requests.get(rawls, headers=head)
        ## TODO add an error msg if we get a 400 etc
        response_dict = json.loads(response.text)
        # then extract the googleProject info
        # google_project_id=response_dict['workspace']['googleProject']
        workspace_name=response_dict['workspace']['name']
        name_output.write(f'{workspace_name}\n')
        workspace_namespace=response_dict['workspace']['namespace']
        namespace_output.write(f'{workspace_namespace}\n')
        return (workspace_namespace, workspace_name)
        #workspace_name=workspace.get_workspace().get('workspace').get('namespace')




def get_sample_sets(workspace_namespace, workspace_name):
    response = requests.get('https://rawls.dsde-prod.broadinstitute.org/api/workspaces/{workspace_namespace}/{workspace_name}/entities?useCache=true')
    sample_sets = response.sample_set
    return sample_sets

# do we want to start and just check for defaults first?
#def get_column_names(workspace_id, workspace_name):


def get_column_values(workspace_id):
    # We need to identify 3 things
    # 1. Sample id field
    # 2. vcf column name
    # 3. vcf index column name

    # We only sample a certain number of rows for each column
    numSamples = 50
    columnSamples = {}

    table_name = "sample" ## TODO is this an assumption that we are making? If the name of the table is sample, then wont the id be sample_id unless override?
    # list_rows is a generator and it uses paging behind the scenes to make
    # this call much more efficient
    for row in table.list_rows(table_name):
        # the id field is special, so handle it like it is...
        if f"{table_name}_id" in columnSamples:
            existingList = columnSamples[f"{table_name}_id"]
            existingList.append(row.name)
        else:
            newList = [row.name];
            columnSamples[f"{table_name}_id"] = newList

        # handle the more general attributes
        for columnName in row.attributes:
            if columnName in columnSamples:
                existingList = columnSamples[columnName]
                existingList.append(row.name)
            else:
                columnSamples[columnName] = [row.attributes[columnName]]

        # done iterating columns
        numSamples -= 1
        if numSamples == 0:
            break


    # time to start gathering some potential rows.
    # how many samples did we actually take?
    numSampledRows = 50 - numSamples
    cutoffPoint = numSampledRows * 0.95
    print(f"Sampled {numSampledRows} rows total. Throwing away any under {cutoffPoint}")


    # match column names that end with "vcf"
    ends_in_vcf_pattern = "^.*vcf$"
    ends_in_vcf = set()

    # match column names that end with "vcf_index"
    ends_in_vcf_index_pattern = "^.*vcf_index$"
    ends_in_vcf_index = set()

    # match column names that contain the word "reblocked"
    contains_reblocked_pattern = ".*reblocked.*"
    contains_reblocked = set()

    # match path that end with ".vcf.gz"
    path_ends_in_vcf_gz_pattern = "^.*\.vcf\.gz$"
    path_ends_in_vcf_gz = set()

    # match path that end with ".vcf.gz.tbi"
    path_ends_in_vcf_gz_tbi_pattern = "^.*\.vcf\.gz\.tbi$"
    path_ends_in_vcf_gz_tbi = set()


    # start sorting the columns that we've seen into buckets to be compared later
    for key in columnSamples:
        print(f"Found key: {key} with {len(columnSamples[key])} entries")
        samplesData = columnSamples[key]
        if len(samplesData) < cutoffPoint:
            # ignoring this completely
            continue

        # ends in vcf?
        result = re.search(ends_in_vcf_pattern, key)
        if result:
            ends_in_vcf.add(key)

        # ends in vcf_index?
        result = re.search(ends_in_vcf_index_pattern, key)
        if result:
            ends_in_vcf_index.add(key)

        # contains the word reblocked?
        result = re.search(contains_reblocked_pattern, key)
        if result:
            contains_reblocked.add(key)

        # has a path that ends in vcf.gz
        result = re.search(path_ends_in_vcf_gz_pattern, samplesData[0])
        if result:
            path_ends_in_vcf_gz.add(key)

        # has a path that ends in vcf.gz.tbi
        result = re.search(path_ends_in_vcf_gz_tbi_pattern, samplesData[0])
        if result:
            path_ends_in_vcf_gz_tbi.add(key)



    print(f"ends_in_vcf: {ends_in_vcf}")
    print(f"ends_in_vcf_index: {ends_in_vcf_index}")
    print(f"contains_reblocked: {contains_reblocked}")

    print(f"path_ends_in_vcf_gz: {path_ends_in_vcf_gz}")
    print(f"path_ends_in_vcf_gz_tbi: {path_ends_in_vcf_gz_tbi}")



    # super simple heuristic: Is there a single entry that ends in vcf
    #for col in ends_in_vcf:


    # and has an analogue


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--workspace_id', type=str,
                        help='The ID of your workspace that holds your sample data',
                        required=True)

    parser.add_argument('--workspace_name', type=str,
                        help='The name of your workspace that holds your sample data',
                        required=False)

    parser.add_argument('--workspace_name_output', type=str,
                        help='The location to write the workspace name to',
                        required=False)

    parser.add_argument('--workspace_namespace_output', type=str,
                        help='The location to write the workspace namespace to',
                        required=False)

    parser.add_argument('--terra_project_id_output', type=str,
                        help='The location to write the terra project id to',
                        required=False)

    parser.add_argument('--attempts_between_pauses', type=int,
                        help='The number of rows in the db that are processed before we pause', default=500)

    args = parser.parse_args()

    # allow this to be overridden, but default it to 500
    if "attempts_between_pauses" in args:
        attempts_between_pauses = args.attempts_between_pauses


    #workspace_name = get_workspace_name(args.workspace_id, args.workspace_name_output, args.workspace_namespace_output)
    # column_names = get_column_values(args.workspace_id, workspace_name)
    column_names = get_column_values(args.workspace_id)

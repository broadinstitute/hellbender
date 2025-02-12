import filecmp
import json
import logging
import os
import re
import pandas as pd
from ast import literal_eval as make_tuple
from typing import List, Optional, Tuple, Set, Dict
import collections

import numpy as np
import pymc as pm

from . import io_consts
from .._version import __version__ as gcnvkernel_version
from ..models.fancy_model import GeneralizedContinuousModel

_logger = logging.getLogger(__name__)


# originally in pymc3.blocking in PyMC3 3.5
VarMap = collections.namedtuple('VarMap', 'var, slc, shp, dtyp')


def read_csv(input_file: str,
             dtypes_dict: Dict[str, object] = None,
             mandatory_columns_set: Set[str] = None,
             comment=io_consts.default_comment_char,
             delimiter=io_consts.default_delimiter_char) -> pd.DataFrame:
    """Opens a file and seeks to the first line that does not start with the comment character,
       checks for mandatory columns in this column-header line, and returns a pandas DataFrame.
       Prefer using this method rather than pandas read_csv, because the comment character will only have an effect 
       when it is present at the beginning of a line at the beginning of the file (pandas can otherwise strip 
       characters that follow a comment character that appears in the middle of a line, which can corrupt sample names).
       Dtypes for columns can be provided, but those for column names that are not known ahead of time will be inferred.

    Args:
        input_file: input file
        dtypes_dict: dictionary of column headers to dtypes; keys will be taken as mandatory columns unless
                     mandatory_columns_set is also provided
        mandatory_columns_set: set of mandatory header columns; must be subset of dtypes_dict keys if provided
        comment: comment character
        delimiter: delimiter character

    Returns:
        pandas DataFrame
    """
    with open(input_file, 'r') as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if not line.startswith(comment):
                fh.seek(pos)
                break
        input_pd = pd.read_csv(fh, delimiter=delimiter, dtype=dtypes_dict)  # dtypes_dict keys may not be present
    found_columns_set = {str(column) for column in input_pd.columns.values}
    assert dtypes_dict is not None or mandatory_columns_set is None, \
        "Cannot specify mandatory_columns_set if dtypes_dict is not specified."
    if dtypes_dict is not None:
        dtypes_dict_keys_set = set(dtypes_dict.keys())
        if mandatory_columns_set is None:
            assert_mandatory_columns(dtypes_dict_keys_set, found_columns_set, input_file)
        else:
            assert mandatory_columns_set.issubset(dtypes_dict_keys_set), \
                "The mandatory_columns_set must be a subset of the dtypes_dict keys."
            assert_mandatory_columns(mandatory_columns_set, found_columns_set, input_file)

    return input_pd


def extract_sample_name_from_header(input_file: str,
                                    comment=io_consts.default_comment_char,
                                    sample_name_header_regexp: str = io_consts.sample_name_header_regexp) -> str:
    """Extracts sample name from header (all lines up to the first line that does not start with the comment character).

    Args:
        input_file: any readable text file
        comment: comment character
        sample_name_header_regexp: the regular expression for identifying the header line that contains
            the sample name

    Returns:
        Sample name
    """
    with open(input_file, 'r') as f:
        while True:
            line = f.readline()
            if not line.startswith(comment):
                break
            match = re.search(sample_name_header_regexp, line, re.M)
            if match is None:
                continue
            groups = match.groups()
            return groups[0]
    raise Exception("Sample name could not be found in \"{0}\"".format(input_file))


def get_sample_name_from_txt_file(input_path: str) -> str:
    """Extract sample name from a text file.

    Args:
        input_path: a path containing the sample name .txt file

    Returns:
        Sample name
    """
    sample_name_file = os.path.join(input_path, io_consts.default_sample_name_txt_filename)
    assert os.path.exists(sample_name_file), \
        "Sample name .txt file could not be found in \"{0}\"".format(input_path)
    with open(sample_name_file, 'r') as f:
        for line in f:
            return line.strip()


def write_sample_name_to_txt_file(output_path: str, sample_name: str):
    """Writes sample name to a text file."""
    with open(os.path.join(output_path, io_consts.default_sample_name_txt_filename), 'w') as f:
        f.write(sample_name + '\n')


def assert_output_path_writable(output_path: str,
                                try_creating_output_path: bool = True):
    """Assert an output path is either writable or can be created upon request.

    Args:
        output_path: the tentative output path
        try_creating_output_path: whether or not try creating the path recursively
            if it does not already exist

    Raises:
        IOError: if the output path is not writable, is not a directory, or does
            not exist and can not be created

    Returns:
        None
    """
    if os.path.exists(output_path):
        if not os.path.isdir(output_path):
            raise IOError("The provided output path \"{0}\" is not a directory".format(output_path))
    elif try_creating_output_path:
        try:
            os.makedirs(output_path)
        except IOError:
            raise IOError("The provided output path \"{0}\" does not exist and can not be created".format(output_path))
    tmp_prefix = "write_tester"
    count = 0
    filename = os.path.join(output_path, tmp_prefix)
    while os.path.exists(filename):
        filename = "{}.{}".format(os.path.join(output_path, tmp_prefix), count)
        count = count + 1
    try:
        filehandle = open(filename, 'w')
        filehandle.close()
        os.remove(filename)
    except IOError:
        raise IOError("The output path \"{0}\" is not writeable".format(output_path))


def write_ndarray_to_tsv(output_file: str,
                         array: np.ndarray,
                         comment=io_consts.default_comment_char,
                         delimiter=io_consts.default_delimiter_char,
                         extra_comment_lines: Optional[List[str]] = None,
                         column_name_str: Optional[str] = None,
                         write_shape_info: bool = True) -> None:
    """Write a vector or matrix ndarray to .tsv file.

    Note:
        Shape and dtype information are stored in the header.

    Args:
        output_file: output .tsv file
        array: array to write to .tsv
        comment: comment character
        delimiter: delimiter character
        extra_comment_lines: (optional) list of extra comment lines to add to the header
        column_name_str: header line (e.g. for representing the ndarray as a table with named columns)
        write_shape_info: if True, ndarray shape info will be written to the header

    Returns:
        None
    """
    array = np.asarray(array)
    assert array.ndim <= 2
    shape = array.shape
    dtype = array.dtype
    header = ""
    if write_shape_info:
        header += compose_sam_comment(io_consts.shape_key_value, repr(shape)) + '\n'
        header += compose_sam_comment(io_consts.type_key_value, str(dtype)) + '\n'
    if extra_comment_lines is not None:
        header += '\n'.join(comment + comment_line for comment_line in extra_comment_lines) + '\n'

    if column_name_str is None:
        header_length = array.shape[1] if array.ndim == 2 else 1
        column_name_str = delimiter.join([io_consts.output_column_prefix + str(i) for i in range(header_length)])
    header += column_name_str + '\n'
    df = pd.DataFrame(array)
    with open(output_file, 'w') as f:
        f.write(header)
        df.to_csv(path_or_buf=f, index=False, header=False, sep=delimiter)


def compose_sam_comment(key: str, value: str) -> str:
    """Compose a SAM style comment string that encodes a key-value pair
    Args:
        key: key string
        value: value string

    Returns:
        A SAM style comment representing the key-value pair

    """
    comment_char = io_consts.default_comment_char
    delim = io_consts.default_delimiter_char
    sep = io_consts.default_key_value_sep
    return comment_char + io_consts.sam_comment_tag + delim + key + sep + value


def parse_sam_comment(comment_line: str) -> Tuple:
    """Parse a SAM style comment

    Args:
        comment_line: a comment string

    Returns:
        Key-value pair represented by a SAM style comment
    """
    match = re.search(io_consts.sam_comment_key_value_regexp, comment_line, re.M)
    if match is None or len(match.groups()) != 2:
        return None, None
    result = match.groups()
    return result[0], result[1]


def read_ndarray_from_tsv(input_file: str,
                          comment=io_consts.default_comment_char,
                          delimiter=io_consts.default_delimiter_char) -> np.ndarray:
    """Reads a vector or matrix ndarray from .tsv file.

    Args:
        input_file: input .tsv file
        comment: comment character
        delimiter: delimiter character

    Returns:
        ndarray
    """
    dtype = None
    shape = None

    with open(input_file, 'r') as f:
        for line in f:
            stripped_line = line.strip()
            if len(stripped_line) == 0:
                continue
            elif stripped_line[0] == comment:
                key, value = parse_sam_comment(stripped_line)
                if key == io_consts.type_key_value:
                    assert dtype is None, "Multiple dtype lines are present in the header of " \
                                          "\"{0}\"".format(input_file)
                    dtype = value
                if key == io_consts.shape_key_value:
                    assert shape is None, "Multiple shape lines are present in the header of " \
                                          "\"{0}\"".format(input_file)
                    shape = make_tuple(value)

    assert dtype is not None and shape is not None, \
        "Shape and dtype information could not be found in the header of " \
        "\"{0}\"".format(input_file)

    df = pd.read_csv(filepath_or_buffer=input_file, sep=delimiter, dtype=dtype, comment=comment)
    return df.values.reshape(shape)


def get_var_map_list_from_mean_field_approx(approx: pm.MeanField) -> List[VarMap]:
    """Extracts the variable-to-linear-array of a PyMC mean-field approximation.

    Args:
        approx: an instance of PyMC mean-field approximation

    Returns:
        A list of VarMap
    """
    # Originally, with PyMC3 3.5, this simply returned a List[pymc3.blocking.VarMap]:
    #     return approx.bij.ordering.vmap
    # However, changes were made to the API and the VarMap class was obviated by the use of Xarray, see:
    #  https://discourse.pymc.io/t/how-to-get-named-means-and-sds-from-advi-fit/11073
    # Unfortunately, this new functionality appears to be somewhat brittle and yields an error in our use case.
    # We instead bring the old VarMap class into this module and recreate the old functionality to
    # preserve our preexisting interfaces.
    var_map_list = []
    for var, slc, shp, dtyp in approx.ordering.values():
        var_map_list.append(VarMap(var, slc, shp, dtyp))
    return var_map_list



def extract_mean_field_posterior_parameters(approx: pm.MeanField) \
        -> Tuple[Set[str], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Extracts mean-field posterior parameters in the right shape and dtype from an instance
    of PyMC mean-field approximation.

    Args:
        approx: an instance of PyMC mean-field approximation

    Returns:
        A tuple (set of variable names,
        map from variable names to their respective Gaussian means,
        map from variable names to their respective Gaussian standard deviations)
    """
    mu_flat_view = approx.mean.get_value()
    std_flat_view = approx.std.eval()
    mu_map = dict()
    std_map = dict()
    var_set = set()
    for vmap in get_var_map_list_from_mean_field_approx(approx):
        var_set.add(vmap.var)
        mu_map[vmap.var] = mu_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
        std_map[vmap.var] = std_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
    return var_set, mu_map, std_map


def write_dict_to_json_file(output_file: str,
                            dict_to_write: Dict,
                            ignored_keys: Set):
    """Writes a dictionary to JSON file.

    Args:
        output_file: output .json file
        dict_to_write: dictionary to write to file
        ignored_keys: a set of keys to ignore
    """
    filtered_dict = {k: v for k, v in sorted(dict_to_write.items()) if k not in ignored_keys}
    with open(output_file, 'w') as fp:
        json.dump(filtered_dict, fp, indent=1)
        fp.write('\n')  # json.dump does not add newline to end of file


def check_gcnvkernel_version_from_json_file(gcnvkernel_version_json_file: str):
    """Reads gcnvkernel version from a JSON file and issues a warning if it is created with a different
    version of the module.

    Args:
        gcnvkernel_version_json_file: input .json file containing gcnvkernel version
    """
    with open(gcnvkernel_version_json_file, 'r') as fp:
        loaded_gcnvkernel_version = json.load(fp)['version']
        if loaded_gcnvkernel_version != gcnvkernel_version:
            _logger.warning("The saved model is created with a different version of gcnvkernel (saved: {0}, "
                            "current: {1}). Backwards compatibility is not guaranteed. Proceed at your own "
                            "risk".format(loaded_gcnvkernel_version, gcnvkernel_version))


def check_gcnvkernel_version_from_path(input_path: str):
    """Reads gcnvkernel version from a path that contains `io_consts.default_gcnvkernel_version_json_filename`,
    reads the gcnvkernel version and issues a warning if it is created with a different version of the module.

    Args:
        input_path:
    """
    check_gcnvkernel_version_from_json_file(
        os.path.join(input_path, io_consts.default_gcnvkernel_version_json_filename))


def _get_mu_tsv_filename(path: str, var_name: str):
    return os.path.join(path, "mu_" + var_name + ".tsv")


def _get_std_tsv_filename(path: str, var_name: str):
    return os.path.join(path, "std_" + var_name + ".tsv")


def _get_singleton_slice_along_axis(array: np.ndarray, axis: int, index: int):
    slc = [slice(None)] * array.ndim
    slc[axis] = index
    return tuple(slc)


def write_mean_field_sample_specific_params(sample_index: int,
                                            sample_posterior_path: str,
                                            approx_var_name_set: Set[str],
                                            approx_mu_map: Dict[str, np.ndarray],
                                            approx_std_map: Dict[str, np.ndarray],
                                            model: GeneralizedContinuousModel,
                                            extra_comment_lines: Optional[List[str]] = None):
    """Writes sample-specific parameters contained in an instance of PyMC mean-field approximation
    to disk.

    Args:
        sample_index: sample integer index
        sample_posterior_path: output path (must be writable)
        approx_var_name_set: set of all variable names in the model
        approx_mu_map: a map from variable names to their respective Gaussian means
        approx_std_map: a map from variable names to their respective Gaussian standard deviations
        model: the generalized model corresponding to the provided mean-field approximation
        extra_comment_lines: (optional) additional comment lines to write to the header of each output file
    """
    sample_specific_var_registry = model.sample_specific_var_registry
    for var_name, var_sample_axis in sample_specific_var_registry.items():
        assert var_name in approx_var_name_set, "A model variable named \"{0}\" could not be found in the " \
                                                "mean-field posterior while trying to write sample-specific " \
                                                "variables to disk".format(var_name)
        mu_all = approx_mu_map[var_name]
        std_all = approx_std_map[var_name]
        mu_slice = np.atleast_1d(mu_all[_get_singleton_slice_along_axis(mu_all, var_sample_axis, sample_index)])
        std_slice = np.atleast_1d(std_all[_get_singleton_slice_along_axis(mu_all, var_sample_axis, sample_index)])

        mu_out_file_name = _get_mu_tsv_filename(sample_posterior_path, var_name)
        write_ndarray_to_tsv(mu_out_file_name, mu_slice, extra_comment_lines=extra_comment_lines)

        std_out_file_name = _get_std_tsv_filename(sample_posterior_path, var_name)
        write_ndarray_to_tsv(std_out_file_name, std_slice, extra_comment_lines=extra_comment_lines)


def write_mean_field_global_params(output_path: str,
                                   approx: pm.MeanField,
                                   model: GeneralizedContinuousModel):
    """Writes global parameters contained in an instance of PyMC mean-field approximation to disk.

    Args:
        output_path: output path (must be writable)
        approx: an instance of PyMC mean-field approximation
        model: the generalized model corresponding to the provided mean-field approximation
    """
    # parse mean-field posterior parameters
    approx_var_set, approx_mu_map, approx_std_map = extract_mean_field_posterior_parameters(approx)

    for var_name in model.global_var_registry:
        assert var_name in approx_var_set, "A model variable named \"{0}\" could not be found in the " \
                                           "mean-field posterior while trying to write global variables " \
                                           "to disk".format(var_name)
        _logger.info("Writing {0}...".format(var_name))
        var_mu = approx_mu_map[var_name]
        var_mu_out_path = _get_mu_tsv_filename(output_path, var_name)
        write_ndarray_to_tsv(var_mu_out_path, var_mu)

        var_std = approx_std_map[var_name]
        var_std_out_path = _get_std_tsv_filename(output_path, var_name)
        write_ndarray_to_tsv(var_std_out_path, var_std)


def read_mean_field_global_params(input_model_path: str,
                                  approx: pm.MeanField,
                                  model: GeneralizedContinuousModel) -> None:
    """Reads global parameters of a given model from saved mean-field posteriors and injects them
    into a provided mean-field instance.

    Args:
        input_model_path: input model path
        approx: an instance of PyMC mean-field approximation to be updated
        model: the generalized model corresponding to the provided mean-field approximation and the saved
            instance
    """
    vmap_list = get_var_map_list_from_mean_field_approx(approx)

    def _update_param_inplace(param, slc, dtype, new_value):
        param[slc] = new_value.astype(dtype).flatten()
        return param

    model_mu = approx.params[0]
    model_rho = approx.params[1]

    for var_name in model.global_var_registry:
        var_mu_input_file = _get_mu_tsv_filename(input_model_path, var_name)
        var_std_input_file = _get_std_tsv_filename(input_model_path, var_name)
        assert os.path.exists(var_mu_input_file) and os.path.exists(var_std_input_file), \
            "Model parameter values for \"{0}\" could not be found in the saved model path while trying " \
            "to read global mean-field parameters".format(var_name)
        _logger.info("Reading model parameter values for \"{0}\"...".format(var_name))
        var_mu = read_ndarray_from_tsv(var_mu_input_file)
        var_std = read_ndarray_from_tsv(var_std_input_file)

        # convert std to rho, see pymc.dist_math.sd2rho
        var_rho = np.log(np.exp(var_std) - 1)
        del var_std

        for vmap in vmap_list:
            if vmap.var == var_name:
                assert var_mu.shape == vmap.shp, \
                    "Loaded mean for \"{0}\" has an unexpected shape; loaded: {1}, " \
                    "expected: {2}".format(var_name, var_mu.shape, vmap.shp)
                assert var_rho.shape == vmap.shp, \
                    "Loaded standard deviation for \"{0}\" has an unexpected shape; loaded: {1}, " \
                    "expected: {2}".format(var_name, var_mu.shape, vmap.shp)
                model_mu.set_value(_update_param_inplace(
                    model_mu.get_value(borrow=True), vmap.slc, vmap.dtyp, var_mu), borrow=True)
                model_rho.set_value(_update_param_inplace(
                    model_rho.get_value(borrow=True), vmap.slc, vmap.dtyp, var_rho), borrow=True)


def read_mean_field_sample_specific_params(input_sample_calls_path: str,
                                           sample_index: int,
                                           sample_name: str,
                                           approx: pm.MeanField,
                                           model: GeneralizedContinuousModel):
    """Reads sample-specific parameters of a given sample from saved mean-field posteriors and injects them
    into a provided mean-field instance.

    Args:
        input_sample_calls_path: path to saved sample-specific posteriors
        sample_index: index of the sample in the current instance of model/approximation
        sample_name: name of the sample in the current instance of model/approximation
            (used to check whether `input_sample_calls_path` actually corresponds to the sample)
        approx: an instance of PyMC mean-field approximation corresponding to the provided model
        model: the generalized model corresponding to the provided mean-field approximation

    Returns:
        None
    """
    path_sample_name = get_sample_name_from_txt_file(input_sample_calls_path)
    assert path_sample_name == sample_name, \
        "The sample name in \"{0}\" does not match the sample name at index {1}; " \
        "found: {2}, expected: {3}. Make sure that the saved posteriors and the current " \
        "task correspond to the same datasets and with the same order/name of samples.".format(
            input_sample_calls_path, sample_index, path_sample_name, sample_name)

    vmap_list = get_var_map_list_from_mean_field_approx(approx)

    def _update_param_inplace(_param: np.ndarray,
                              _var_slice: slice,
                              _var_shape: Tuple,
                              _sample_specific_loaded_value: np.ndarray,
                              _var_sample_axis: int,
                              _sample_index: int) -> np.ndarray:
        """Updates the ndarray buffer of the shared parameter tensor according to a given sample-specific
        parameter for a given sample index.

        Args:
            _param: ndarray buffer of the shared parameter tensor (i.e. `mu` or `rho`)
            _var_slice: the slice that of `_param` that yields the full view of the sample-specific
                parameter to be updated
            _var_shape: full shape of the sample-specific parameter to be updated
            _sample_specific_loaded_value: new single-sample slice of the sample-specific parameter
                to be updated
            _var_sample_axis: the sample-index axis in the full view of the sample-specific
                parameter to be updates
            _sample_index: sample index

        Returns:
            updated `_param`
        """
        sample_specific_var = _param[_var_slice].reshape(_var_shape)
        sample_specific_var[_get_singleton_slice_along_axis(
            sample_specific_var, _var_sample_axis, _sample_index)] = _sample_specific_loaded_value[:]
        return _param

    # reference to mean-field posterior mu and rho
    model_mu = approx.params[0]
    model_rho = approx.params[1]

    for var_name, var_sample_axis in model.sample_specific_var_registry.items():
        var_mu_input_file = _get_mu_tsv_filename(input_sample_calls_path, var_name)
        var_std_input_file = _get_std_tsv_filename(input_sample_calls_path, var_name)
        assert os.path.exists(var_mu_input_file) and os.path.exists(var_std_input_file), \
            "Model parameter values for \"{0}\" could not be found in the provided calls " \
            "path \"{1}\"".format(var_name, input_sample_calls_path)

        var_mu = read_ndarray_from_tsv(var_mu_input_file)
        var_std = read_ndarray_from_tsv(var_std_input_file)

        # convert std to rho, see pymc.dist_math.sd2rho
        var_rho = np.log(np.exp(var_std) - 1)
        del var_std

        # update mu and rho
        for vmap in vmap_list:
            if vmap.var == var_name:
                model_mu.set_value(_update_param_inplace(
                    model_mu.get_value(borrow=True), vmap.slc, vmap.shp, var_mu,
                    var_sample_axis, sample_index), borrow=True)
                model_rho.set_value(_update_param_inplace(
                    model_rho.get_value(borrow=True), vmap.slc, vmap.shp, var_rho,
                    var_sample_axis, sample_index), borrow=True)


def write_gcnvkernel_version(output_path: str):
    """Writes the current gcnvkernel version as a JSON file to a given path.

    Args:
        output_path: path to write the gcnvkernel version

    Returns:
        None
    """
    # write gcnvkernel version
    write_dict_to_json_file(
        os.path.join(output_path, io_consts.default_gcnvkernel_version_json_filename),
        {'version': gcnvkernel_version}, set())


def assert_mandatory_columns(mandatory_columns_set: Set[str],
                             found_columns_set: Set[str],
                             input_tsv_file: str):
    """Asserts that a given .tsv file contains a set of mandatory header columns.

    Note:
        The set of header columns found in the .tsv file must be provided. `input_tsv_file` is only used
        for generating exception messages.
    Args:
        mandatory_columns_set: set of mandatory header columns
        found_columns_set: set of header columns found in the .tsv file
        input_tsv_file: path to the .tsv file in question

    Returns:
        None
    """
    not_found_set = mandatory_columns_set.difference(found_columns_set)
    assert len(not_found_set) == 0, "The following mandatory columns could not be found in \"{0}\"; " \
                                    "cannot continue: {1}".format(input_tsv_file, not_found_set)


def assert_files_are_identical(input_file_1: str, input_file_2: str):
    """Asserts that two given files are bit identical."""
    assert os.path.isfile(input_file_1), "Cannot find {0}.".format(input_file_1)
    assert os.path.isfile(input_file_2), "Cannot find {0}.".format(input_file_2)
    assert filecmp.cmp(input_file_1, input_file_2, shallow=False), \
        "The following two files are expected to be identical: {0}, {1}".format(input_file_1, input_file_2)

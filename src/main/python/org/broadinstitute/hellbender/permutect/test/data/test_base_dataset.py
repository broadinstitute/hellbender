import tempfile
import permutect.data.base_dataset as ds
import torch

from permutect import utils


def test_line_to_tensor():
    line1 = "1.0 1.1 1.2 1.3"
    tensor1 = ds.line_to_tensor(line1)

    # allow for tensor rounding error
    assert torch.max(torch.tensor([1.0, 1.1, 1.2, 1.3]) - tensor1).item() < 0.001


def test_read_integers():
    line1 = "1 2 3 4"
    integers = ds.read_integers(line1)

    # allow for tensor rounding error
    assert list(integers) == [1, 2, 3, 4]


def test_read_2d_tensor():
    tmp = tempfile.NamedTemporaryFile()

    with open(tmp.name, 'w') as f:
        lines = [
            "1 2 3 4 5\n",
            "6 7 8 9 10\n",
            "11 12 13 14 15\n",
            "NOTHING\n",
            "10 20 30\n",
            "40 50 60\n"
        ]

        f.writelines(lines)

    with open(tmp.name) as f:
        tensor1 = ds.read_2d_tensor(f, 3)
        f.readline()
        tensor2 = ds.read_2d_tensor(f, 2)

        assert list(tensor1.size()) == [3, 5]
        assert list(tensor2.size()) == [2, 3]

        assert tensor1.tolist() == [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]]


def test_read_data():
    tmp = tempfile.NamedTemporaryFile()

    with open(tmp.name, 'w') as f:
        lines = [
            "UNLABELED\n",
            "1:12807, C->T\n",
            "GGAGAGGCTTCGATGCCCCTC\n",
            "0.192 0.000 0.000 1.000 1.000 1.000 1.000 1.000 1.000\n"
            "5 2 10 2\n"
            "23 30 1 0 42 21 362 42 320 0 0\n",
            "27 32 1 0 9 21 290 9 281 0 0\n",
            "24 30 0 0 15 21 353 15 338 0 0\n",
            "24 30 0 0 15 21 353 15 338 0 0\n",
            "24 30 0 0 15 21 353 15 338 0 0\n",
            "24 32 0 0 23 42 351 77 274 0 0\n",
            "23 31 0 0 35 30 350 65 285 0 0\n",
            "36 4 26 2\n",
            "0.879\n",
            "4.55\n",
            "ARTIFACT\n",
            "1:13079, C->G\n",
            "CCAGCTGGGTCGACAGACAGG\n",
            "0.113 0.045 3.000 0.833 1.000 1.000 1.000 1.000 1.000\n",
            "5 1 10 2\n",
            "27 24 0 1 8 21 290 281 9 0 0\n",
            "27 24 0 1 8 21 290 281 9 0 0\n",
            "27 24 0 1 8 21 290 281 9 0 0\n",
            "27 24 0 1 8 21 290 281 9 0 0\n",
            "27 24 0 1 8 21 290 281 9 0 0\n",
            "23 31 1 1 4 21 346 341 5 0 0\n",
            "78 4 75 2\n",
            "12.5\n",
            "-73.3\n"
        ]

        f.writelines(lines)

    data = list(ds.read_data(tmp.name))
    assert len(data) == 2
    assert data[0].label == utils.Label.UNLABELED
    assert torch.max(data[0].get_info_tensor_1d() - torch.tensor([0.192, 0.000, 0.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000] + [1, 0, 0])).item() < 0.001

    assert data[1].reads_2d.size()[0] == 6

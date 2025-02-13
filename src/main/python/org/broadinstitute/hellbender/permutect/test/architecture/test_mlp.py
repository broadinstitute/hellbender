import torch
from permutect import utils

from permutect.architecture.mlp import MLP


# test with artificial data where a*x = 0 is a perfect linear separator
def test_linearly_separable_data():
    input_dim = 3
    num_samples = 1000
    a = torch.rand(input_dim, 1)
    x = torch.rand(num_samples, input_dim)
    y = (torch.sign(torch.matmul(x,a)) + 1)/2   # labels are 0 / 1

    layer_sizes = [input_dim, 1]
    model = MLP(layer_sizes)

    loss_func = torch.nn.BCEWithLogitsLoss(reduction='mean')
    optimizer = torch.optim.Adam(model.parameters())
    loss_list = []

    num_epochs = 10000
    for epoch in range(num_epochs):
        prediction = model.forward(x)
        loss = loss_func(prediction, y)
        utils.backpropagate(optimizer, loss)
        loss_list.append(loss.item())

    assert loss_list[-1] < 0.01
    assert loss_list[1000] < loss_list[0]
    assert loss_list[2000] < loss_list[1000]
    assert loss_list[3000] < loss_list[2000]

    pred = torch.sign(torch.sigmoid(model.forward(x))-0.5)
    lab = torch.sign(y-0.5)

    errors = torch.sum(torch.abs((pred - lab)/2)).item()

    # we should have perfect accuracy
    assert errors < 1


# test with annular data where y = 1 when  1/3 < norm(x) < 2/3
def test_annular_data():
    input_dim = 3
    num_samples = 1000
    x = torch.randn(num_samples, input_dim)/torch.sqrt(torch.tensor([input_dim]))

    norms = torch.norm(x, dim=1)
    y = (torch.sign(norms - 0.33) * torch.sign(0.66 - norms) + 1)/2

    layer_sizes = [input_dim, 5, 5, 5, 5, 1]
    model = MLP(layer_sizes)

    loss_func = torch.nn.BCEWithLogitsLoss(reduction='mean')
    optimizer = torch.optim.Adam(model.parameters())
    loss_list = []

    num_epochs = 10000
    for epoch in range(num_epochs):
        prediction = model.forward(x)
        loss = loss_func(torch.squeeze(prediction), y)
        utils.backpropagate(optimizer, loss)
        loss_list.append(loss.item())

    assert loss_list[-1] < 0.2

    pred = torch.squeeze(torch.sign(torch.sigmoid(model.forward(x)) - 0.5))
    lab = torch.sign(y - 0.5)

    errors = torch.sum(torch.abs((pred - lab) / 2)).item()

    # we should have perfect accuracy
    assert errors < 100


# A single hidden layer should suffice to perfectly classify a 2D XOR pattern where
# the 1st and 3rd quadrants are labeled 0 and the 2nd and 4th are 1
def test_xor_data():
    num_samples = 1000
    x = torch.randn(num_samples, 2)
    y = (torch.sign(x[:, 0]*x[:, 1]) + 1)/2   # labels are 0 / 1

    # Use a single hidden layer.  Really only 2 neurons are needed, but more speeds convergence
    layer_sizes = [2, 4, 1]
    model = MLP(layer_sizes)

    loss_func = torch.nn.BCEWithLogitsLoss(reduction='mean')
    optimizer = torch.optim.Adam(model.parameters())
    loss_list = []

    num_epochs = 10000
    for epoch in range(num_epochs):
        prediction = model.forward(x)
        loss = loss_func(prediction.squeeze(), y)
        utils.backpropagate(optimizer, loss)
        loss_list.append(loss.item())

    assert loss_list[-1] < 0.15
    assert loss_list[1000] < loss_list[0]
    assert loss_list[2000] < loss_list[1000]
    assert loss_list[3000] < loss_list[2000]

    pred = torch.sign(torch.sigmoid(model.forward(x).detach().squeeze())-0.5)
    lab = torch.sign(y-0.5)

    errors = torch.sum(torch.abs((pred - lab)/2)).item()

    # we should have near-perfect accuracy
    assert errors < num_samples / 50

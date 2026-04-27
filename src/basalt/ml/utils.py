#!/usr/bin/env python

"""
Utility functions for BASALT deep learning components.

This module provides helpers for downloading model weights, normalising
feature matrices, managing logging directories, training and evaluating
PyTorch models, and computing focal loss.
"""

import os
import sys
from glob import glob

from torch import nn
import torch
from torch.nn import functional as F
from tqdm import tqdm
import numpy as np

import requests


def download_model(url, local_dir=None):
    """
    Download a model file from a remote URL into the BASALT weight directory.

    If ``local_dir`` is not provided, the function will use the directory
    specified by the ``BASALT_WEIGHT`` environment variable. The file is
    downloaded with a progress bar and will not be downloaded again if it
    already exists.
    """
    if local_dir is None:
        user_dir = os.path.expanduser('~')
        # Default cache directory could be under the user home, but here
        # we rely on an explicit BASALT_WEIGHT environment variable to
        # avoid implicit side effects.
        local_dir = os.environ.get("BASALT_WEIGHT")
        if not local_dir:
            raise EnvironmentError(
                "BASALT_WEIGHT environment variable is not set. "
                "Please configure a directory for BASALT model weights."
            )
    local_path = f"{local_dir}/{os.path.basename(url)}"

    if os.path.exists(local_path):
        print(f"File already exists in {local_path}.")
        return local_dir

    os.makedirs(local_dir, exist_ok=True)
    print(f"File will be saved in {local_path}.")

    # Stream the file to disk with a progress bar to support large files.
    response = requests.get(url, stream=True)
    response.raise_for_status()
    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024
    progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True)

    with open(local_path, 'wb') as f:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            f.write(data)

    progress_bar.close()
    return local_dir

def del_best_ckpt(log_dir):
    """
    Remove all checkpoint files containing ``'best'`` in the given log directory.

    This is typically used to clean up older best checkpoints when
    re-running experiments.
    """
    for i in glob(os.path.join(log_dir, '*best*')):
        os.remove(i)

# def norm(t, type='mmn'):
#     if type == 'mmn':
#         return (t - np.min(t, axis=0)) / (np.max(t, axis=0) - np.min(t, axis=0))
#     elif type == 'absmmn':
#         t = np.abs(t-np.mean(t, axis=0))
#         return t/np.max(t, axis=0)
#     else:
#         return t

def norm(t, type='mmn'):
    """
    Normalize a NumPy array along columns using one of the supported strategies.

    Parameters
    ----------
    t : np.ndarray
        Input tensor/array to be normalized.
    type : {'mmn', 'absmmn'}, optional
        Normalization type. ``'mmn'`` performs min-max normalization;
        ``'absmmn'`` performs min-max normalization on the absolute
        deviation from the mean. Any other value will return the input
        unchanged.

    Returns
    -------
    np.ndarray
        Normalized array.
    """
    if type == 'mmn':
        return (t - np.min(t, axis=0)) / (np.max(t, axis=0) - np.min(t, axis=0) + 1e-8)
    elif type == 'absmmn':
        t = np.abs(t-np.mean(t, axis=0))
        return t/(np.max(t, axis=0) + 1e-8)
    else:
        return t


def get_log_dir(set='train', comment=''):
    """
    Create and return a unique log directory for the current run.

    Parameters
    ----------
    set : str, optional
        Tag for the run type, e.g. ``'train'`` or ``'val'``.
    comment : str, optional
        Additional string appended to the log directory name.

    Returns
    -------
    str
        Absolute path to the created log directory.

    Raises
    ------
    FileExistsError
        If a directory with the same name already exists.
    """
    log_dir = os.path.join('runs', set, comment)
    if os.path.exists(log_dir):
        raise FileExistsError
    os.makedirs(log_dir, exist_ok=True)
    return log_dir


def train_one_epoch(model, optimizer, data_loader, device, epoch, weight, loss_fun):
    """
    Train the model for a single epoch.

    This function iterates over the provided data loader, performs forward
    and backward passes, and updates model parameters using the specified
    optimizer and loss function.
    """
    model.train()
    if loss_fun == 'ce':
        loss_function = torch.nn.CrossEntropyLoss(weight=torch.tensor([weight, 1-weight], dtype=torch.float).to(device))
        # loss_function = torch.nn.CrossEntropyLoss()
    elif loss_fun == 'fl':
        loss_function = focal_loss(alpha=weight, num_classes=2)
    mean_loss = torch.zeros(1).to(device)
    optimizer.zero_grad()

    data_loader = tqdm(data_loader, file=sys.stdout)

    for step, data in enumerate(data_loader):
        x, labels = data

        pred = model(x.to(device))

        loss = loss_function(pred, labels.to(device))
        loss.backward()
        mean_loss = (mean_loss * step + loss.detach()) / (step + 1)  # update mean losses

        data_loader.desc = "[train epoch {}] mean loss {}".format(epoch, round(mean_loss.item(), 3))

        if not torch.isfinite(loss):
            print('WARNING: non-finite loss, ending training ', loss)
            sys.exit(1)

        optimizer.step()
        optimizer.zero_grad()

    return mean_loss.item()


@torch.no_grad()
def evaluate(model, data_loader, device, epoch, weight, loss_fun):
    """
    Evaluate the model on a validation set.

    Returns per-class accuracies for class 0 and class 1 to better capture
    performance on imbalanced datasets.
    """
    model.eval()

    if loss_fun == 'ce':
        loss_function = torch.nn.CrossEntropyLoss(weight=torch.tensor([weight, 1-weight], dtype=torch.float).to(device))
        # loss_function = torch.nn.CrossEntropyLoss()
    elif loss_fun == 'fl':
        loss_function = focal_loss(alpha=weight, num_classes=2)

    tr, pr = [], []

    mean_loss = torch.zeros(1).to(device)
    data_loader = tqdm(data_loader, file=sys.stdout)

    for step, data in enumerate(data_loader):
        x, labels = data
        pred = model(x.to(device))

        loss = loss_function(pred, labels.to(device))
        mean_loss = (mean_loss * step + loss.detach()) / (step + 1)  # update mean losses
        data_loader.desc = "[val epoch {}] mean loss {}".format(epoch, round(mean_loss.item(), 3))

        pred = torch.max(pred, dim=1)[1]
        # sum_num += torch.eq(pred, labels.to(device)).sum()
        tr += list(labels.detach().cpu().numpy().flatten())
        pr += list(pred.detach().cpu().numpy().flatten())

    tr = np.array(tr)
    pr = np.array(pr)
    # print(accuracy_score(tr, pr))
    mask_0 = tr == 0
    tr_0, pr_0 = tr[mask_0], pr[mask_0]
    acc_0 = np.sum(tr_0==pr_0)/tr_0.shape[0]
    mask_1 = tr == 1
    tr_1, pr_1 = tr[mask_1], pr[mask_1]
    acc_1 = np.sum(tr_1==pr_1)/tr_1.shape[0]

    return acc_0, acc_1

def save_confusion_mat(log_dir, matrix, classes=None, title=None):
    """
    Plot and save a confusion matrix figure to the given log directory.
    """
    from matplotlib import pyplot as plt
    from sklearn.metrics import ConfusionMatrixDisplay
    plt.figure()
    disp = ConfusionMatrixDisplay(confusion_matrix=matrix, display_labels=classes)
    disp.plot()
    if title:
        plt.title(title)
    name = title or 'confusion_mat'
    plt.savefig(f'{log_dir}/{name}.png')


@torch.no_grad()
def evaluate_ensemble(model, data_loader, device):
    """
    Evaluate an ensemble model and return the predicted labels.

    This helper assumes that labels are not required and only predictions
    are needed downstream (e.g. for voting or thresholding).
    """
    model.eval()


    tr, pr = [], []

    for step, data in enumerate(data_loader):
        x = data
        pred = model(x.to(device))
        pred = torch.max(pred, dim=1)[1]
        # sum_num += torch.eq(pred, labels.to(device)).sum()
        # tr += list(labels.detach().cpu().numpy().flatten())
        pr += list(pred.detach().cpu().numpy().flatten())

    # tr = np.array(tr)
    pr = np.array(pr)

    # print(accuracy_score(tr, pr))

    # mask_0 = tr == 0
    # tr_0, pr_0 = tr[mask_0], pr[mask_0]
    # acc_0 = np.sum(tr_0==pr_0)/tr_0.shape[0]
    # mask_1 = tr == 1
    # tr_1, pr_1 = tr[mask_1], pr[mask_1]
    # acc_1 = np.sum(tr_1==pr_1)/tr_1.shape[0]
    #
    # print(f'Contaminated total: {np.sum(mask_1)}\n'
    #       f'Contaminated correctly removed: {np.sum(tr_1==pr_1)}\n'
    #       f'Contaminated wrong remained: {np.sum(mask_1)-np.sum(tr_1==pr_1)}\n'
    #       f'Real total: {np.sum(mask_0)}\n'
    #       f'Real wrong removed: {np.sum(mask_0)-np.sum(tr_0==pr_0)}\n'
    #       f'Real correctly remained: {np.sum(tr_0==pr_0)}\n')

    return pr


class focal_loss(nn.Module):
    """
    Focal loss for binary or multi-class classification.

    This loss is designed to address class imbalance by down-weighting
    well-classified examples and focusing training on hard examples.
    """

    def __init__(self, alpha=0.25, gamma=2, num_classes=3, size_average=True):
        super(focal_loss, self).__init__()
        self.size_average = size_average
        if isinstance(alpha, list):
            assert len(alpha) == num_classes
            self.alpha = torch.Tensor(alpha)
        else:
            assert alpha < 1
            self.alpha = torch.zeros(num_classes)
            self.alpha[0] += alpha
            self.alpha[1:] += (1 - alpha)

        self.gamma = gamma

    def forward(self, preds, labels):
        """
        Compute focal loss given predictions and ground-truth labels.
        """
        preds = preds.view(-1, preds.size(-1))
        self.alpha = self.alpha.to(preds.device)
        preds_logsoft = F.log_softmax(preds, dim=1)  # log_softmax
        preds_softmax = torch.exp(preds_logsoft)  # softmax

        # Implement NLL loss by gathering the probabilities of the true class
        preds_softmax = preds_softmax.gather(1, labels.view(-1, 1))
        preds_logsoft = preds_logsoft.gather(1, labels.view(-1, 1))
        self.alpha = self.alpha.gather(0, labels.view(-1))

        # (1 - p_t) ** gamma term from focal loss
        loss = -torch.mul(torch.pow((1 - preds_softmax), self.gamma), preds_logsoft)

        loss = torch.mul(self.alpha, loss.t())
        if self.size_average:
            loss = loss.mean()
        else:
            loss = loss.sum()
        return loss

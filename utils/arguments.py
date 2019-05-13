"""Define the complete arguments for the program."""

import argparse


def get_args():
    """Return the complete arguments."""
    parser = argparse.ArgumentParser(
        description='Evaluation arguments for LTR elements.')

    # Personal settings
    pred = '/Users/tefa/Documents/Transposones/Data/results/yeast_noSTR.txt'
    gt = '/Users/tefa/Documents/Transposones/Data/yeast/' + \
         'gt_transposons.txt'

    parser.add_argument('--pred', type=str, default=pred,
                        help='path to the prediction file')
    parser.add_argument('--gt', type=str, default=gt,
                        help='path to the groundtruth file')

    return parser.parse_args()

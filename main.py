"""Evaluation code for LTR."""
import numpy as np
import matplotlib.pyplot as plt
from models.transposon import Transposon
from utils.arguments import get_args


def process_file(filename, mode='gt'):
    """Process the groundtruth/predictions file."""
    print('Processing {} as groundtruth'.format(filename))
    f = open(filename, 'r')
    LTRs = {}
    act_name = ''
    for line in f.readlines():
        data = line.split('\t')
        if act_name != data[0]:
            act_LTRs = []
            act_name = data[0]
            act_LTRs.append(Transposon(data[0], data[1], data[2], data[3]))
            LTRs[act_name] = act_LTRs
        else:
            LTRs[act_name].append(Transposon(data[0], data[1], data[2],
                                             data[3]))
    f.close()
    return LTRs


def get_single_instance_results(gts, preds, thresh):
    """Calculate the tp, fp, fn for a given chromosome."""
    tp = 0
    IoUs = []
    gt_index = []
    pred_index = []
    for x, pred in enumerate(preds):
        for y, gt in enumerate(gts):
            act_overlap = pred.get_overlap(gt)
            IoU = act_overlap/(len(pred) + len(gt) - act_overlap)
            if IoU >= thresh:
                IoUs.append(IoU)
                gt_index.append(y)
                pred_index.append(x)

    IoUs = np.argsort(IoUs)[::-1]
    if len(IoUs) == 0:
        # No matches
        return 0, len(preds), len(gts)

    gt_match_idx = []
    pred_match_idx = []
    for idx in IoUs:
        gt_idx = gt_index[idx]
        pr_idx = pred_index[idx]
        if (gt_idx not in gt_match_idx) and (pr_idx not in pred_match_idx):
            gt_match_idx.append(gt_idx)
            pred_match_idx.append(pr_idx)
    tp = len(gt_match_idx)
    fp = len(preds) - len(pred_match_idx)
    fn = len(gts) - len(gt_match_idx)
    return tp, fp, fn


def calculate_metrics(gt, pred):
    """Calculate the cumulative PR values for different thresholds."""
    thresholds = np.linspace(0.01, 1, 100)
    precisions = []
    recalls = []
    names = list(set(gt.keys()).union(set(pred.keys())))
    for thresh in thresholds:
        cum_tp, cum_fp, cum_fn = 0, 0, 0
        for seq_name in names:
            if seq_name not in pred:
                cum_fn += len(gt[seq_name])
                continue
            elif seq_name not in gt:
                cum_fp += len(pred[seq_name])
                continue
            act_gt = gt[seq_name]
            act_pred = pred[seq_name]
            tp, fp, fn = get_single_instance_results(act_gt, act_pred, thresh)
            cum_tp += tp
            cum_fp += fp
            cum_fn += fn
            print(cum_tp, cum_fp, cum_fn)
        print('-----------')
        print('threshold', thresh)
        print('precision', cum_tp/(cum_tp + cum_fp + 1e-9))
        print('recall', cum_tp/(cum_tp + cum_fn + 1e-9))
        print('-----------')
        precisions.append(cum_tp/(cum_tp + cum_fp + 1e-9))
        recalls.append(cum_tp/(cum_tp + cum_fn + 1e-9))
    return precisions, recalls


def plot_PR(precisions, recalls):
    """Plot precision-recall curve."""
    print(precisions)
    print(recalls)
    plt.plot(recalls, precisions, 'ro')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision Recall Curve for Transposon Detection')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()


def main():
    """Run main function for the program."""
    args = get_args()
    args_dict = vars(args)
    print(args_dict)
    print('Argument list to program')
    print('\n'.join(['--{0} {1}'.format(arg, args_dict[arg])
                    for arg in args_dict]))
    print('\n')

    gt = process_file(args.gt)
    pred = process_file(args.pred, mode='pred')

    precisions, recalls = calculate_metrics(gt, pred)

    plot_PR(precisions, recalls)


# Run main
if __name__ == '__main__':
    main()

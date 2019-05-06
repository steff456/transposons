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
        if mode == 'gt':
            data = line.split('\t')
            if act_name != data[0]:
                act_LTRs = []
                act_name = data[0]
                act_LTRs.append(Transposon(data[0], data[1], data[2]))
                LTRs[act_name] = act_LTRs
            else:
                LTRs[act_name].append(Transposon(data[0], data[1], data[2]))
        else:
            data = line.split('\t')
            if act_name != data[0]:
                act_LTRs = []
                act_name = data[0]
                act_LTRs.append(Transposon(data[0], data[1], data[2]))
                LTRs[act_name] = act_LTRs
            else:
                LTRs[act_name].append(Transposon(data[0], data[1], data[2]))
    f.close()
    return LTRs


def get_single_instance_results(gts, preds, thresh):
    """Calculate the tp, fp, fn for a given chromosome."""
    tp = 0
    for pred in preds:
        for gt in gts:
            act_overlap = pred.get_overlap(gt)
            if act_overlap >= thresh:
                tp += 1
    fp = len(preds) - tp
    fn = len(gts) - tp
    return tp, fp, fn


def calculate_metrics(gt, pred):
    """Calculate the cumulative PR values for different thresholds."""
    thresholds = np.linspace(0.5, 1, 20)
    precisions = []
    recalls = []
    for thresh in thresholds:
        cum_tp, cum_fp, cum_fn = 0, 0, 0
        for seq_name in gt.keys():
            act_gt = gt[seq_name]
            act_pred = pred[seq_name]
            tp, fp, fn = get_single_instance_results(act_gt, act_pred, thresh)
            cum_tp += tp
            cum_fp += fp
            cum_fn += fn
        precisions.append(cum_tp/(cum_tp + cum_fp + 1e-9))
        recalls.append(cum_tp/(cum_tp + cum_fn + 1e-9))
    return precisions, recalls


def plot_PR(precisions, recalls):
    """Plot precision-recall curve."""
    plt.plot(recalls, precisions)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision Recall Curve for Transposon Detection')
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

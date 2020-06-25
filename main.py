"""Evaluation code for LTR."""
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
from heapq import heappush, heappop

from utils.arguments import get_args
from models.transposon import Transposon


def process_file(filename, mode='gt'):
    """Process the groundtruth/predictions file."""
    print('Processing {} as {}'.format(filename, mode))
    f = open(filename, 'r')
    LTRs = {}
    act_name = ''
    for line in f.readlines():
        data = line.split('\t')
        if act_name != data[0]:
            act_LTRs = []
            act_name = data[0]
            if mode == 'gt':
                act_LTRs.append(Transposon(data[0], data[1], data[2], 0))
            else:
                act_LTRs.append(Transposon(data[0], data[1], data[2], data[3]))
            LTRs[act_name] = act_LTRs
        else:
            if mode == 'gt':
                act_LTRs.append(Transposon(data[0], data[1], data[2], 0))
            else:
                act_LTRs.append(Transposon(data[0], data[1], data[2], data[3]))
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


def calculate_metrics(gt, pred, thresh=0.5, verbose=True):
    """Calculate the cumulative PR values for different thresholds."""
    names = list(set(gt.keys()).union(set(pred.keys())))
    cum_tp, cum_fp, cum_fn, total_p, total_gt = 0, 0, 0, 0, 0
    for seq_name in names:
        if seq_name not in pred:
            cum_fn += len(gt[seq_name])
            total_gt += len(gt[seq_name])
            continue
        elif seq_name not in gt:
            cum_fp += len(pred[seq_name])
            total_p += len(pred[seq_name])
            continue
        act_gt = gt[seq_name]
        act_pred = pred[seq_name]
        tp, fp, fn = get_single_instance_results(act_gt, act_pred, thresh)
        if False:
            precision = calculate_precision(tp, fp)
            recall = calculate_recall(tp, fn)
            print('----------- {} -----------'.format(seq_name))
            print('TP:', tp, 'FP:', fp, 'FN:', fn)
            print('precision', precision)
            print('recall', recall)
            print('F-measure', calculate_fmeasure(precision, recall))
            print('-----------')
        cum_tp += tp
        cum_fp += fp
        cum_fn += fn
        total_p += len(pred[seq_name])
        total_gt += len(gt[seq_name])
    precision = calculate_precision(cum_tp, cum_fp)
    recall = calculate_recall(cum_tp, cum_fn)
    fmeasure = calculate_fmeasure(precision, recall)
    if verbose:
        print('-----------')
        print('TP:', cum_tp, 'FP:', cum_fp, 'FN:', cum_fn)
        print('Total pred:', total_p, 'Total gt:', total_gt)
        print('threshold', thresh)
        print('precision', precision)
        print('recall', recall)
        print('F-measure', fmeasure)
        print('-----------')
    return precision, recall, fmeasure, cum_fp, cum_tp


def calculate_precision(tp, fp):
    """Calculate precision metric given TP and FP."""
    return tp/(tp + fp + 1e-9)


def calculate_recall(tp, fn):
    """Calculate recall metric given TP and FN."""
    return tp/(tp + fn + 1e-9)


def calculate_fmeasure(precision, recall):
    """Calculate F-measure."""
    return (2*precision*recall)/(precision + recall + 1e-9)


def get_preds_scores_map(preds):
    """Get the scores for the predictions in ascending order."""
    genome_scores = {}
    for name in preds:
        scores = []
        for act in preds[name]:
            heappush(scores, (act.score, act))
        genome_scores[name] = scores
    return genome_scores


def calculate_different_recalls_single_thresh(gt, preds, thresh, scores):
    """Calculate precision at different recalls for the predictions."""
    precisions = []
    recalls = []
    false_positives = []
    true_positives = []
    for score in scores:
        # pred = get_preds_scores_map(preds)
        pred = preds
        total_p = []
        total_r = []
        for name in pred:
            for transposon in pred[name]:
                if transposon.score < score:
                    heappop(pred[name])
        precision, recall, fm, fp, tp = calculate_metrics(
            gt, pred, thresh=thresh)
        total_p.append(precision)
        total_r.append(recall)
        print('-------- Score {} ---------'.format(score))
        print('Precision: {}'.format(precision))
        print('Recall: {}'.format(recall))
        print('Fmeasure: {}'.format(fm))
        precisions.append(precision)
        recalls.append(recall)
        false_positives.append(fp)
        true_positives.append(tp)
    return precisions, recalls, false_positives, true_positives


def plot_PR(precisions, recalls, threshs=0.5):
    """Plot precision-recall curve."""
    # for precision, recall, thresh in zip(precisions, recalls, threshs):
    plt.plot(recalls, precisions, label=threshs)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision Recall Curve for Transposon Detection')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(loc='best')
    plt.show()


def plot_ROC(fpr, tpr):
    """Plot ROC curve."""
    roc_auc = metrics.auc(fpr, tpr)
    display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc,
                                      estimator_name='TransposonFinder')
    display.plot()
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
    # threshs = [0.5, 0.6, 0.7, 0.8, 0.9]
    scores = [0, 2000, 3000, 5000, 8000, 9000, 10000]
    scores = np.linspace(0, 16000, 5)
    total_p, total_r, fps, tps = (
        calculate_different_recalls_single_thresh(gt, pred, 0.5, scores))
    plot_PR(total_p, total_r, [0.5])
    plot_ROC(fps, tps)


# Run main
if __name__ == '__main__':
    main()

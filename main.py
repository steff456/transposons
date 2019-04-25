"""Evaluation code for LTR."""
from models.transposon import Transposon
from utils.arguments import get_args
from statistics import mean 


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


def compare(gt, pred):
    """Find overlap in the gt and predicted sites."""
    precision = []  # Avg precision in all the genome
    recall = []  # Avg recall in all the genome
    total_gt = 0  # Total number of transposons in the gt file

    for seq_name in gt.keys():
        act_gt = gt[seq_name]
        if seq_name not in pred:
            total_gt += len(act_gt)
            continue
        act_pred = pred[seq_name]
        gt_index = 0
        pred_index = 0
        while gt_index < len(act_gt) - 1:
            if pred_index > len(act_pred) - 1:
                break
            gt_TE = act_gt[gt_index]
            pred_TE = act_pred[pred_index]
            # Calculate all the metrics for the act pair if there's an overlap
            overlap = gt_TE.is_overlap(pred_TE)
            if overlap:
                size_o = gt_TE.get_overlap(pred_TE)
                recall.append(size_o/len(gt_TE))
                precision.append(size_o/len(pred_TE))
                pred_index += 1
                gt_index += 1
            # Advance condition
            elif pred_TE.first < gt_TE.first:
                pred_index += 1
            else:
                gt_index += 1
        total_gt += len(act_gt)

    avg_recall = mean(recall)
    avg_precision = mean(precision)
    f_measure = (2*avg_precision*avg_recall)/(avg_precision+avg_recall)

    print('Recall: {}'.format(avg_recall))
    print('Precision: {}'.format(avg_precision))
    print('F-measure: {}'.format(f_measure))


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

    compare(gt, pred)


# Run main
if __name__ == '__main__':
    main()

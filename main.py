"""Evaluation code for LTR."""
from models.transposon import Transposon
from utils.arguments import get_args
import pdb


def process_file(filename, mode='gt'):
    """Process the groundtruth/predictions file."""
    print('Processing {} as groundtruth'.format(filename))
    f = open(filename, 'r')
    LTRs = {}
    act_name = ''
    for line in f.readlines():
        if mode == 'gt':
            if 'long_terminal_repeat' in line or 'LTR' in line:
                data = line.split('\t')
                if act_name != data[0]:
                    act_LTRs = []
                    act_name = data[0]
                    act_LTRs.append(Transposon(data[0], data[1], data[2]))
                    LTRs[act_name] = act_LTRs
                else:
                    LTRs[act_name].append(Transposon(data[0],
                                                     data[1], data[2]))
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
    found = 0  # Total predictions that are in the gt
    total_gt = 0  # Total TEs expected by the gt
    total_pred = 0  # Total TEs predicted

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
            if abs(gt_TE.first - pred_TE.first) < 200:
                found += 1
                pred_index += 1
                gt_index += 1
                print('Difference of size {} vs {}'.format(len(gt_TE),
                                                           len(pred_TE)))
            # Advance condition
            elif pred_TE.first < gt_TE.first:
                pred_index += 1
            else:
                gt_index += 1
        total_gt += len(act_gt)
        total_pred += len(act_pred)

    print('Found: {} \nTotal Gt: {} \nTotal Pred: {}'.format(found,
                                                             total_gt,
                                                             total_pred))


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

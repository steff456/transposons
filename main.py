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
                    act_LTRs[act_name].append(Transposon(data[0],
                                                         data[1], data[2]))
        else:
            data = line.split('\t')
            LTRs.append(Transposon(data[0], data[1], data[2]))
    f.close()
    return LTRs


def compare(gt, pred):
    """Find overlap in the gt and predicted sites."""
    for seq_name in gt.keys():
        act_gt = gt[seq_name]
        act_pred = pred[seq_name]
        pdb.set_trace()


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

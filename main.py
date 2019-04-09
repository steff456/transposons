"""Evaluation code for LTR."""
from models.transposon import Transposon
from utils.arguments import get_args


def process_file(filename, mode='gt'):
    """Process the groundtruth/predictions file."""
    print('Processing {} as groundtruth'.format(filename))
    f = open(filename, 'r')
    LTRs = []
    for line in f.readlines():
        if mode == 'gt':
            if 'long_terminal_repeat' in line or 'LTR' in line:
                data = line.split('\t')
                LTRs.append(Transposon(data[0], data[1], data[2]))
        else:
            data = line.split('\t')
            LTRs.append(Transposon(data[0], data[1], data[2]))
    f.close()
    return LTRs


def compare(gt, pred):
    """Find overlap in the gt and predicted sites."""


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

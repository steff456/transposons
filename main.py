"""Evaluation code for LTR."""
from utils.arguments import get_args


def main():
    """Run main function for the program."""
    args = get_args()
    args_dict = vars(args)
    print(args_dict)
    print('Argument list to program')
    print('\n'.join(['--{0} {1}'.format(arg, args_dict[arg])
                    for arg in args_dict]))
    print('\n')


# Run main
if __name__ == '__main__':
    main()

#Create in file given parameters

import argparse
import os

def main():

    parser = argparse.ArgumentParser(description='Write active noise input file with given input parameters.')
    parser.add_argument('indir',
                        help='Directory where input file will be written.')
    parser.add_argument('--dim', default='2')
    parser.add_argument('--dx', default='1.0')
    parser.add_argument('--nx', default='16')
    parser.add_argument('--ny', default='16')
    #parser.add_argument('--nz', default='16')
    parser.add_argument('--tau', default='1.0')
    parser.add_argument('--Lambda', default='1.0')
    parser.add_argument('--D', default='1.0')
    parser.add_argument('--do_fft', default='1')

    parser.add_argument('--dt', default='1.0')

    parser.add_argument('--output_dir', default='../active-noise-results')
    parser.add_argument('--freq', default='1')

    args = parser.parse_args()
    print(args.__dict__)

    try:
        os.makedirs(args.indir)
        print('Made directory.')
    except OSError as e:
        print('Directory exists')

    print('Writing input file...')
    with open(args.indir + '/active_noise_nx=%d_ny=%d_tau=%f_lambda=%f.in' % (int(args.nx),int(args.ny),float(args.tau),float(args.Lambda)), 'w') as f:
        f.write('#Generator\n')
        f.write('dim = %s\n' % args.dim)
        f.write('dx = %s\n' % args.dx)
        f.write('nx = %s\n' % args.nx)
        f.write('ny = %s\n' % args.ny)
        #f.write('nz = %s\n' % args.nz)
        f.write('tau = %s\n' % args.tau)
        f.write('lambda = %s\n' % args.Lambda)
        f.write('D = %s\n' % args.D)
        f.write('do_fft = %s\n' % args.do_fft)
        f.write('\n')
        f.write('#Solver\n')
        f.write('dt = %s\n' % args.dt)
        f.write('\n')
        f.write('#Observer\n')
        f.write('output_dir = %s\n' % args.output_dir)
        f.write('freq = %s\n' % args.freq)

if __name__=='__main__':
    main()

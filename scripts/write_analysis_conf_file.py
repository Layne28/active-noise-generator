#Create analysis conf file given parameters

import argparse
import os

def main():

    parser = argparse.ArgumentParser(description='Write active noise conf file with given input parameters.')
    parser.add_argument('confdir',
                        help='Directory where conf file will be written.')
    parser.add_argument('--output_dir', default='../active-noise-results/')
    parser.add_argument('--input_dir', default='../active-noise-results/')
    parser.add_argument('--nx', default='16')
    parser.add_argument('--ny', default='16')
    parser.add_argument('--nz', default='16')
    parser.add_argument('--tau', default='0.01')
    parser.add_argument('--Lambda', default='8.0')
    parser.add_argument('--nframes', default='10000')
    parser.add_argument('--nchunks', default='10')
    parser.add_argument('--frame_interval', default='10')
    parser.add_argument('--rmax', default='8.0')
    parser.add_argument('--delta_x', default='1.0')

    parser.add_argument('--tmax', default='10.0')

    args = parser.parse_args()
    print(args.__dict__)

    try:
        os.makedirs(args.confdir)
        print('Made directory.')
    except OSError as e:
        print('Directory exists')

    print('Writing conf file...')
    with open(args.confdir + '/noise_analysis_nx=%d_ny=%d_nz=%d_tau=%f_lambda=%f.conf' % (int(args.nx),int(args.ny),int(args.nz),float(args.tau),float(args.Lambda)), 'w') as f:
        f.write('output_dir = %s\n' % args.output_dir)
        f.write('input_dir = %s\n' % args.input_dir)
        f.write('nx = %s\n' % args.nx)
        f.write('ny = %s\n' % args.ny)
        f.write('nz = %s\n' % args.nz)
        f.write('tau = %s\n' % args.tau)
        f.write('lambda = %s\n' % args.Lambda)
        f.write('nframes = %s\n' % args.nframes)
        f.write('nchunks = %s\n' % args.nchunks)
        f.write('frame_interval = %s\n' % args.frame_interval)

if __name__=='__main__':
    main()
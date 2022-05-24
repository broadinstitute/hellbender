import argparse
import sys


def scale_xy_bed_values(input_file, output_file, x_scale_factor, y_scale_factor):
    input_bed = open(input_file, 'r')
    output_bed = open(output_file, 'w')

    fail = False

    if x_scale_factor < 1.0:
        print(f"Error: illegal X chromosome weight scale factor {x_scale_factor}; scale factor value must be >= 1.0")
        fail = True

    if y_scale_factor < 1.0:
        print(f"Error: illegal Y chromosome weight scale factor {y_scale_factor}; scale factor value must be >= 1.0")
        fail = True

    if fail:
        sys.exit(1)

    while True:
        line = input_bed.readline()
        if not line:
            break

        line = line.rstrip('\n')

        if line.startswith('chrX') or line.startswith('chrY'):
            scale_factor = x_scale_factor if line.startswith('chrX') else y_scale_factor
            fields = line.split('\t')
            weight = fields[-1]
            fields[-1] = str(int(int(weight) * scale_factor))
            output_bed.write('\t'.join(fields) + '\n')
        else:
            output_bed.write(line + '\n')

    output_bed.close()


def parse_args():
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Scale X and Y BED values for more uniform extract shard runtimes')
    parser.add_argument('--input', type=str, help='Input BED file', required=True)
    parser.add_argument('--output', type=str, help='Output BED file', required=True)
    parser.add_argument('--xscale', type=float, help='X chromosome scaling factor', required=True)
    parser.add_argument('--yscale', type=float, help='Y chromosome scaling factor', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    scale_xy_bed_values(args.input, args.output, args.xscale, args.yscale)

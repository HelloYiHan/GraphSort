import os
#os.chdir('GraphSort') 先把要计算的文件都放到这个程序所在的目录下
import argparse
import torch
import graphsort_core #把原来的graphsort.py改名为graphsort_core.py
import read_dataset
#import graphsort_test
import os.path as osp
import torch_geometric.data
import subprocess


#parameters
parser = argparse.ArgumentParser(description='Arguments of GraphSort')
parser.add_argument('--input', '-i', type=str, help = 'Expression data of samples', required = True)
parser.add_argument('--type', '-t', type=str, help = 'Data type, rnaseq OR microarray)',required = True)
parser.add_argument('--output', '-o', type=str, help = 'File name of output', default = 'graphsort_out.txt')
parser.add_argument('--batch_size', '-b', type=int, help = 'Batch size for computation', default = '10')
parser.add_argument('--device', '-d', type=str, help = 'Computation device, gpu OR cpu', default = 'cpu')
args = parser.parse_args()

print(args.input)
print(args.type)
print(args.output)
print(args.batch_size)
print(args.device)

#select device
if args.device == 'cpu':
    device = torch.device('cpu')
elif args.device == 'gpu' and torch.cuda.is_available():
    device = torch.device('cuda')
elif args.device == 'gpu' and not torch.cuda.is_available():
    raise RuntimeError('GPU was selected but not available in the machine.')
else:
    raise RuntimeError('Device argument (--device or -d) must be gpu or cpu.')

path = os.path.dirname(os.path.abspath(__file__))

#preprocess
if args.type == 'rnaseq':
        train_file = path + '/rem_bat_eff_dat_n1000.txt'
elif args.type == 'microarray':
        train_file = path + '/rem_bat_eff_dat_microarray.txt'
else:
    raise RuntimeError('Data type argument (--type or -t) must be rnaseq or microarray.')

subprocess.call(['chmod', 'u+x', path + '/preprocessing.R'])
subprocess.call([path + '/preprocessing.R', args.input, train_file, args.type])

#load model
model = graphsort_core.Net(num_features = 1, y_size = 8).to(device)
if args.type == 'rnaseq':
    if device.type == "cuda":
        model.load_state_dict(torch.load(path + '/trained_models/state1627195299007epoch60.pkl',map_location='cuda'))
    else:
        model.load_state_dict(torch.load(path + '/trained_models/state1627195299007epoch60.pkl',map_location='cpu'))
elif args.type =='microarray':
    if device.type == "cuda":
        model.load_state_dict(torch.load('./trained_models/.pkl',map_location='cuda'))
    else:
        model.load_state_dict(torch.load('./trained_models/.pkl',map_location='cpu'))

#estimation
dataset_input = read_dataset.two_dim_data(zippath= 'InputFile_' + args.input + '.zip', root = osp.join('.', 'data', 'InputFile_' + args.input), name = 'InputFile_' + args.input, use_node_attr = True)
loader_input = torch_geometric.data.DataLoader(dataset_input, batch_size=args.batch_size)
graphsort_core.estimate(model, loader_input, device, args.output)
print("GraphSort estimation done")
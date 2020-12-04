import os
import argparse
import torch
import graphsort_core
import read_dataset
import os.path as osp
import torch_geometric.data
import subprocess
import numpy as np

#1.arguments
parser = argparse.ArgumentParser(description='Arguments of GraphSort')
parser.add_argument('--input', '-i', type=str, help = 'Expression data of samples', required = True)
parser.add_argument('--type', '-t', type=str, help = 'Data type, rnaseq, microarray, OR pancreatic)',required = True)
parser.add_argument('--output', '-o', type=str, help = 'File name of output', default = 'graphsort_out.txt')
parser.add_argument('--batch_size', '-b', type=int, help = 'Batch size for computation', default = 10)
parser.add_argument('--device', '-d', type=str, help = 'Computation device, gpu OR cpu', default = 'cpu')
parser.add_argument('--size', '-s', type=str, help = 'Size of preprocessing data, 1k OR 2k', default = '1k')
args = parser.parse_args()

#2.select device
if args.device == 'cpu':
    device = torch.device('cpu')
elif args.device == 'gpu' and torch.cuda.is_available():
    device = torch.device('cuda')
elif args.device == 'gpu' and not torch.cuda.is_available():
    raise RuntimeError('GPU was selected but not available in the machine.')
else:
    raise RuntimeError('Device argument (--device or -d) must be gpu or cpu.')

path = os.path.dirname(os.path.abspath(__file__))

#3.preprocess
print('Start preprocessing')
if args.type == 'rnaseq':
    if args.size == '1k':
        train_file = path + '/rem_bat_eff_dat_n1000.txt'
    elif args.size == '2k':
        train_file = path + '/rem_bat_eff_dat_n2000.txt'
    else:
        raise RuntimeError('Size of preprocessing data (--size or -s) must be 1k or 2k.')
elif args.type == 'microarray':
    train_file = path + '/rem_bat_eff_dat_microarray.txt'
elif args.type == 'pancreatic':
    train_file = path + '/rem_bat_eff_dat_pancreatic.txt'
else:
    raise RuntimeError('Data type argument (--type or -t) must be rnaseq or microarray.')

subprocess.call(['chmod', 'u+x', path + '/preprocessing.R'])
subprocess.call([path + '/preprocessing.R', args.input, train_file, args.type])
print('Preprocessing done')

#4.load model
model = graphsort_core.Net(num_features = 1, y_size = 8).to(device)
if args.type == 'rnaseq':
    if device.type == "cuda":
        model.load_state_dict(torch.load(path + '/trained_models/graphsort_rnaseq_model.pkl',map_location='cuda'))
    else:
        model.load_state_dict(torch.load(path + '/trained_models/graphsort_rnaseq_model.pkl',map_location='cpu'))
elif args.type =='microarray':
    if device.type == "cuda":
        model.load_state_dict(torch.load(path + '/trained_models/graphsort_microarray_model.pkl',map_location='cuda'))
    else:
        model.load_state_dict(torch.load(path + '/trained_models/graphsort_microarray_model.pkl',map_location='cpu'))
elif args.type == 'pancreatic':
    if device.type == "cuda":
        model.load_state_dict(torch.load(path + '/trained_models/graphsort_pancreatic_model.pkl',map_location='cuda'))
    else:
        model.load_state_dict(torch.load(path + '/trained_models/graphsort_pancreatic_model.pkl',map_location='cpu'))

#5.estimation
print('Start estimation')
dataset_input = read_dataset.two_dim_data(zippath= 'InputFile_' + args.input + '.zip', root = osp.join('.', 'data', 'InputFile_' + args.input), name = 'InputFile_' + args.input, use_node_attr = True)
loader_input = torch_geometric.data.DataLoader(dataset_input, batch_size=args.batch_size)
if args.type == 'rnaseq':
    graphsort_core.estimateRNASEQ(model, loader_input, device, args.output, args.input)
elif args.type == 'microarray':
    graphsort_core.estimateMicroarray(model, loader_input, device, args.output, args.input)
elif args.type == 'pancreatic':
    graphsort_core.estimatePancreatic(model, loader_input, device, args.output, args.input)
print("GraphSort estimation done")
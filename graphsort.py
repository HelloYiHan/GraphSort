import numpy as np
import pandas as pd
import os.path as osp
import shutil
import torch
import torch.nn.functional as F
from torch_geometric.nn import ARMAConv, TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp

class Net(torch.nn.Module):
    def __init__(self, y_size, num_features=1):
        self.num_features = num_features
        self.y_size  = y_size
        super(Net, self).__init__()
        
        self.conv1 = ARMAConv(self.num_features, 128)
        self.bn1 = torch.nn.BatchNorm1d(128)
        self.pool1 = TopKPooling(128, ratio=0.8)
        
        self.conv2 = ARMAConv(128, 128)
        self.bn2 = torch.nn.BatchNorm1d(128)
        self.pool2 = TopKPooling(128, ratio=0.8)
        
        self.conv3 = ARMAConv(128, 128)
        self.bn3 = torch.nn.BatchNorm1d(128)
        self.pool3 = TopKPooling(128, ratio=0.8)
        
        self.conv4 = ARMAConv(128,128)
        self.bn4 = torch.nn.BatchNorm1d(128)
        self.pool4 = TopKPooling(128, ratio=0.8)
        
        self.conv5 = ARMAConv(128,128)
        self.bn5 = torch.nn.BatchNorm1d(128)
        self.pool5 = TopKPooling(128, ratio=0.8)
        
        self.conv6 = ARMAConv(128,128)
        self.bn6 = torch.nn.BatchNorm1d(128)
        self.pool6 = TopKPooling(128, ratio=0.8)
        
        self.conv7 = ARMAConv(128, 128)
        self.bn7 = torch.nn.BatchNorm1d(128)
        self.pool7 = TopKPooling(128, ratio=0.8)
        
        self.conv8 = ARMAConv(128,128)
        self.bn8 = torch.nn.BatchNorm1d(128)
        self.pool8 = TopKPooling(128, ratio=0.8)
        
        self.conv9 = ARMAConv(128,128)
        self.bn9 = torch.nn.BatchNorm1d(128)
        self.pool9 = TopKPooling(128, ratio=0.8)
        
        self.conv10 = ARMAConv(128,128)
        self.bn10 = torch.nn.BatchNorm1d(128)
        self.pool10 = TopKPooling(128, ratio=0.8)
        
        self.lin1 = torch.nn.Linear(2560,1280)
        self.lin2 = torch.nn.Linear(1280, self.y_size)
        
        self.act = torch.nn.PReLU()
        
    def min_max_norm(self,A):
        A -= A.min(1, keepdim=True)[0]
        A /= A.max(1, keepdim=True)[0]
        A /= torch.sum(A, dim=1,keepdim=True)
        return(A)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        
        x = self.conv1(x, edge_index)
        x = self.bn1(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        
        x = self.conv2(x, edge_index)
        x = self.bn2(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool2(x, edge_index, None, batch)
        x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        
        x = self.conv3(x, edge_index)
        x = self.bn3(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool3(x, edge_index, None, batch)
        x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
    
        x = self.conv4(x, edge_index)
        x = self.bn4(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool4(x, edge_index, None, batch)
        x4 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = self.conv5(x, edge_index)
        x = self.bn5(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool5(x, edge_index, None, batch)
        x5 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = self.conv6(x, edge_index)
        x = self.bn6(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool6(x, edge_index, None, batch)
        x6 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        
        x = self.conv7(x, edge_index)
        x = self.bn7(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool7(x, edge_index, None, batch)
        x7 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        
        x = self.conv8(x, edge_index)
        x = self.bn8(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool8(x, edge_index, None, batch)
        x8 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = self.conv9(x, edge_index)
        x = self.bn9(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool9(x, edge_index, None, batch)
        x9 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = self.conv10(x, edge_index)
        x = self.bn10(x)
        x = self.act(x)
        x, edge_index, _, batch, _, _ = self.pool10(x, edge_index, None, batch)
        x10 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
               
        x = torch.cat([x1, x2, x3, x4, x5, x6, x7, x8, x9,x10],dim=1)
        x = self.act(self.lin1(x))
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.act(self.lin2(x))
        x = self.min_max_norm(x)
        
        return x



def test_65133(loader):
    model.eval()
    cor_all = 0
    for data in loader:
        data = data.to(device)
        x = model(data)
        
        xb = x[:,0] + x[:,1]
        xtcd4 = x[:,3]
        xtcd8 = x[:,4]
        xmon = x[:,5]
        xnk = x[:,7]
        x5 = torch.cat((xb.unsqueeze(1),xtcd4.unsqueeze(1),xtcd8.unsqueeze(1),xmon.unsqueeze(1),xnk.unsqueeze(1)),1)
            
        yb = data.y[:,0] + data.y[:,1]
        ytcd4 = data.y[:,3] + data.y[:,4] + data.y[:,5]
        ytcd8 = data.y[:,2]
        ymon = data.y[:,8]
        ynk = data.y[:,7]
        y5 = torch.cat((yb.unsqueeze(1),ytcd4.unsqueeze(1),ytcd8.unsqueeze(1),ymon.unsqueeze(1),ynk.unsqueeze(1)),1)
        
        
        cor_sample = pearson(x5, y5)
        cor_all += torch.sum(cor_sample,dim=0).item()
        
    cor_all /= len(loader.dataset)

    return cor_all


def test_106898(loader):
    model.eval()
    cor_all = 0
    for data in loader:
        data = data.to(device)
        x = model(data)
        
        xb = x[:,0] + x[:,1] + x[:,2]
        xtcd4 = x[:,3]
        xtcd8 = x[:,4]
        xmon = x[:,5]
        xdc = x[:,6]
        xnk = x[:,7]
        x6 = torch.cat((xb.unsqueeze(1),xtcd4.unsqueeze(1),xtcd8.unsqueeze(1),xmon.unsqueeze(1),xdc.unsqueeze(1),xnk.unsqueeze(1)),1)
            
        yb = data.y[:,0] + data.y[:,1] + data.y[:,2] + data.y[:,3] + data.y[:,4]
        ytcd4 = data.y[:,5]
        ytcd8 = data.y[:,6]
        ymon = data.y[:,10] + data.y[:,11] + data.y[:,12]
        ydc = data.y[:,8] + data.y[:,9]
        ynk = data.y[:,13]
        y6 = torch.cat((yb.unsqueeze(1),ytcd4.unsqueeze(1),ytcd8.unsqueeze(1),ymon.unsqueeze(1),ydc.unsqueeze(1),ynk.unsqueeze(1)),1)
        
        cor_sample = pearson(x6, y6)
        cor_all += torch.sum(cor_sample,dim=0).item()

    cor_all /= len(loader.dataset)

    return cor_all
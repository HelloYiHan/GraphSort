import torch

def test_rna_seq(model, loader, device):
    model.eval()
    cor_all = 0
    for data in loader:
        data = data.to(device)
        x = model(data)
        cor = pearson(x[:,:7], data.y[:,:7])
        cor_all += torch.sum(cor,dim=0).item()
    return cor_all / len(loader.dataset)


def test_gse65133(model, loader, device):
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


def test_gse106898(model, loader, device):
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


def pearson(x, y):
    x = x - torch.mean(x, dim=1)[:,None]
    y = y - torch.mean(y, dim=1)[:,None]
    cor = torch.sum(x * y,dim=1) / (torch.sqrt(torch.sum(x ** 2,dim=1)) * torch.sqrt(torch.sum(y ** 2,dim=1)))
    return cor


def estimate(loader, rFile):
    model.eval()
    for data in loader:
        data = data.to(device)
        x = model(data)
        with open(rFile, "a") as fr:
            np.savetxt(fr,x.cpu().detach().numpy())
    return 0
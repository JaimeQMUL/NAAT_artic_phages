import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from src.seq2topt.functions import *
import timeit
import torch
from src.seq2topt.model import MultiAttModel
import esm
from src.biotools.fasta_tools import *
from src.seq2topt import *

# LOAD tOPT MODEL
# path to model weights
topt_pth = '../src/seq2topt/model_topt_window.3_r2.0.57.pth';# you can download the weights from release v1.0.0
# define device
device = torch.device('cpu')
# define seq2Topt model
emb_dim= 320; window=3; n_head = 4; n_RD = 4;
model = MultiAttModel( emb_dim, window, n_head, n_RD)
model.to(device);
# load weights into the model
model.load_state_dict(torch.load( topt_pth, map_location=device  ));
model.eval();

# LOAD ESM-2
# We use esm2_t6_8M_UR50D, 6 layers, dimension size = 320.
esm2_model, alphabet = esm.pretrained.esm2_t6_8M_UR50D() # 6 layers
esm2_model = esm2_model.to(device)
esm2_batch_converter = alphabet.get_batch_converter()

def predict(seq):
    # seq: input protein sequence
    inputs = [('Temp', seq)]
    batch_labels, batch_strs, batch_tokens = esm2_batch_converter(inputs)
    batch_tokens = batch_tokens.to(device=device, non_blocking=True)
    with torch.no_grad():
        emb = esm2_model(batch_tokens, repr_layers=[6], return_contacts=False)
    emb = emb["representations"][6]
    emb = emb.transpose(1,2)
    emb = emb.to(device)
    with torch.no_grad():
        preds = model( emb )
    return preds

Topt_max=120

# Get top hits accessions
top_hits_accs=pd.read_csv('top_hits/metadata/top_hits.csv')

tOPT=[]
for acc in top_hits_accs['accession']:
    seq=ExtractSequence(acc, 'top_hits/sequences/top_hits.fasta')
    temp=int(predict(seq)*Topt_max)
    tOPT.append(temp)



# Load TM model
# path to model weights
topt_pth = '../src/seq2topt/model_tm_window.3_r2.0.76.pth' # you can download the weights from release v1.0.0
# define device
device = torch.device('cpu')
# define seq2Topt model
emb_dim= 320; window=3; n_head = 4; n_RD = 4;
model = MultiAttModel( emb_dim, window, n_head, n_RD)
model.to(device);
# load weights into the model
model.load_state_dict(torch.load( topt_pth, map_location=device  ));
model.eval();


# We use esm2_t6_8M_UR50D, 6 layers, dimension size = 320.
esm2_model, alphabet = esm.pretrained.esm2_t6_8M_UR50D() # 6 layers
esm2_model = esm2_model.to(device)
esm2_batch_converter = alphabet.get_batch_converter()

# We use esm2_t6_8M_UR50D, 6 layers, dimension size = 320.
esm2_model, alphabet = esm.pretrained.esm2_t6_8M_UR50D() # 6 layers
esm2_model = esm2_model.to(device)
esm2_batch_converter = alphabet.get_batch_converter()

TM=[]
for acc in top_hits_accs['accession']:
    seq=ExtractSequence(acc, 'top_hits/sequences/top_hits.fasta')
    temp=int(predict(seq)*Topt_max)
    TM.append(temp)

df = pd.DataFrame({
    'accession': top_hits_accs['accession'],
    'tOPT': tOPT,
    'TM': TM
})

df.to_csv('top_hits/topt/top_hits_tp.csv', index=False)
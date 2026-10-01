#!/usr/bin/env python3
"""Deterministic figure implementation. Frozen input values; no analyses."""
from pathlib import Path
import os, sys, hashlib, json, shutil, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from PIL import Image, ImageDraw

BASE=Path(os.environ.get('NEUROTRACE_PROJECT_ROOT', '.')).resolve()
P=BASE/'figure_source_data'
OUT=Path(os.environ.get('NEUROTRACE_FIGURE_OUTPUT', str(BASE/'generated_figures')))
MODS=['NTM1_ASD_up','NTM2_ASD_down','NTM3_ASD_signed']
MC=['#7B3294','#2166AC','#1B7837']
MODULE_LABELS=['NTM1 upregulated','NTM2 downregulated','NTM3 mixed']
def module_label(p): return MODULE_LABELS[int(mod(p)[-1])-1]
CC={'positive':'#7B3294','negative':'#2166AC'}
STAGES=['early_prenatal','mid_prenatal','late_prenatal','early_postnatal','childhood_adolescence','adulthood']
SL=['EPr','MPr','LPr','EPost','C/A','Adult']
SHORT=['EPr','MPr','LPr','EPost','C/A','Adult']
STAGE_KEY='EPr, early prenatal; MPr, mid prenatal; LPr, late prenatal; EPost, early postnatal; C/A, childhood–adolescence; Adult, adulthood.'
INK='#27313A'; GRAY='#7D858C'; LIGHT='#E8ECEE'
CMAP=LinearSegmentedColormap.from_list('teal', ['#F4F7F6','#A3C8C5','#17636A'])
DIV=LinearSegmentedColormap.from_list('signed',['#356A9B','#F8F8F6','#B85D42'])
CMAP.set_bad('#E4E6E8'); DIV.set_bad('#E4E6E8')
METHODS=['DIRECT_NATIVE_PROJECTION','DUAL_HEAD_GENE_FIRST_PPR','GENE_ONLY_PPR_STAGE3','GRAPH_OT','HETERO_PPR_STANDARD','LEGACY_NEUROTRACE_SIM_PROXY','MEAN_SIGNATURE','RANDOM','SIGNED_GENE_FIRST_PPR_SINGLE_HEAD','SIMPLE_OT','STAGE_CORRELATION','UNSIGNED_GENE_PPR']
ML=['Native projection','Dual-head PPR','Gene-only PPR','Graph OT','Heterogeneous PPR','Simulation proxy','Mean signature','Random','Single-head signed PPR','Simple OT','Stage correlation','Unsigned PPR']
CODES=['NP','DH','GP','GO','HP','SP','MS','RD','SH','SO','SC','UP']
COND={'baseline':'Baseline','high_batch':'High batch','composition_shift':'Composition shift','high_dropout':'High dropout','false_prior':'False prior','combined_stress':'Combined stress','weak_signal':'Weak signal','low_overlap':'Low overlap','composition_confounded':'Composition confound','false_prior_stress':'False-prior stress','combined_hard':'Combined challenge'}
READS=[]; EXPORTS=[]; TEXTS=[]; FIG=''
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7.5,'legend.fontsize':7.5,'axes.linewidth':.65,'axes.edgecolor':GRAY,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white','lines.linewidth':1.3,'lines.markersize':4,'axes.spines.top':False,'axes.spines.right':False})

def load(fn, block=''):
    d=pd.read_csv(P/fn,sep='\t')
    READS.append(dict(figure=FIG,block=block,input_file=str(P/fn),rows=len(d),columns=';'.join(d.columns),sha256=hashlib.sha256((P/fn).read_bytes()).hexdigest()))
    return d

def fig(name, h=8.0):
    global FIG
    FIG=name
    return plt.figure(figsize=(8.3,h))

def ax(f,box): return f.add_axes(box)
def title(f, letter, name, x, y):
    f.text(x,y,letter,weight='bold',fontsize=11,va='bottom')
    f.text(x+.031,y,name,weight='bold',fontsize=9,va='bottom')
def note(f,text,x=.08,y=.025,size=7): f.text(x,y,text,fontsize=size,va='bottom')
def clean(a,grid='y'):
    a.set_axisbelow(True)
    a.grid(axis=grid,color=LIGHT,lw=.55)
    a.tick_params(length=3,width=.6,pad=3)
def zero(a,vertical=False):
    (a.axvline if vertical else a.axhline)(0,color='#9A9EA1',lw=.7,zorder=0)
def mod(p): return p.split('_')[0].split(' ')[0]
def color(p): return MC[int(mod(p)[-1])-1]
def fmt(v, n=2): return '—' if pd.isna(v) else f'{v:.{n}f}'
def fdr(v): return f'{v:.5g}'
def threshold_legend(f,y=.96):
    f.legend([Line2D([],[],color=INK,ls='--',marker='o',mfc='white'),Line2D([],[],color=INK,ls='-',marker='s')],['Top 200','Top 500'],loc='center right',bbox_to_anchor=(.96,y),ncol=2,frameon=False,handlelength=2.3,columnspacing=1.5)
def module_legend(f,y=.025):
    f.legend([Line2D([],[],color=c,marker='o') for c in MC],['NTM1','NTM2','NTM3'],loc='center',bbox_to_anchor=(.5,y),ncol=3,frameon=False)
def heat(a,values,rows,cols,lo=None,hi=None,cmap=CMAP,annot=True,fs=7):
    v=np.array(values,dtype=float)
    im=a.pcolormesh(np.arange(v.shape[1]+1),np.arange(v.shape[0]+1),np.ma.masked_invalid(v),cmap=cmap,vmin=lo,vmax=hi,edgecolors='white',linewidth=.65)
    a.set_xlim(0,v.shape[1]); a.set_ylim(v.shape[0],0)
    a.set_xticks(np.arange(v.shape[1])+.5,cols)
    a.set_yticks(np.arange(v.shape[0])+.5,rows)
    a.tick_params(length=0,pad=5)
    for s in a.spines.values(): s.set_visible(False)
    if annot:
        norm=im.norm
        for i in range(v.shape[0]):
            for j in range(v.shape[1]):
                val=v[i,j]; rgba=im.cmap(norm(val)) if np.isfinite(val) else (1,1,1,1)
                lum=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
                a.text(j+.5,i+.5,fmt(val),ha='center',va='center',fontsize=fs,color='white' if lum<.5 else INK)
    return im
def cb(f,im,rect,label,vertical=False):
    c=f.colorbar(im,cax=ax(f,rect),orientation='vertical' if vertical else 'horizontal')
    c.outline.set_visible(False); c.ax.tick_params(length=2,labelsize=7,pad=2); c.set_label(label,size=7,labelpad=3)
    if hasattr(c,'solids') and c.solids is not None: c.solids.set_rasterized(False)

def profile(a,d,col='magnitude_cosine',hue='program',thresholds=(200,500),ylim=(-.32,.34),labels=True):
    for key,g in d.groupby(hue,sort=False):
        c=CC[key] if hue=='channel' else color(key)
        for t in thresholds:
            q=g[g.top_n==t].sort_values('stage_order')
            if q.empty or q[col].isna().all(): continue
            ls,marker=('--','o') if t==200 else ('-','s')
            if len(thresholds)==4:
                ls={50:':',100:'-.',200:'--',500:'-'}[t]; marker={50:'v',100:'^',200:'o',500:'s'}[t]
            a.plot(range(6),q[col],color=c,ls=ls,marker=marker,mfc='white' if t!=500 else c,ms=3.4,lw=1.25)
    zero(a); a.axvline(2.5,color=LIGHT,lw=1)
    a.set_xticks(range(6),SHORT); a.set_xlim(-.25,5.25); a.set_ylim(*ylim)
    if labels: a.set_ylabel('Cosine similarity')
    clean(a)
    if d[col].isna().all(): a.text(.5,.62,'Channel absent',ha='center',va='center',transform=a.transAxes,color=GRAY,fontsize=8)


def channel_composition(a,d):
    rows=[(p,t) for p in MODS for t in [200,500]]
    for i,(p,t) in enumerate(rows):
        left=0
        for ch in ['positive','negative']:
            q=d[(d.program==p)&(d.top_n==t)&(d.channel==ch)].iloc[0]; share=q.channel_mass_share_of_total_abs_weight
            a.barh(i,share,left=left,color=CC[ch],height=.58,edgecolor='white',lw=.6)
            if share>0: a.text(left+share/2,i,f'{share:.0%}',ha='center',va='center',fontsize=7,color='white')
            left+=share
        pos=d[(d.program==p)&(d.top_n==t)&(d.channel=='positive')].iloc[0].n_restart_genes
        neg=d[(d.program==p)&(d.top_n==t)&(d.channel=='negative')].iloc[0].n_restart_genes
        a.text(1.04,i,f'{int(pos)} / {int(neg)}',va='center',fontsize=7)
    a.set_yticks(range(6),[f'{mod(p)} · Top {t}' for p,t in rows]); a.set_ylim(5.6,-.6); a.set_xlim(0,1.33)
    a.set_xticks([0,.5,1],['0','0.5','1']); a.set_xlabel('Share of absolute restart weight')
    a.text(1.035,-.72,'Genes + / −',fontsize=7); clean(a,'x')

def s1():
    name='Supplementary_Figure_S1'; f=fig(name,8.3)
    d=load('Supplementary_FigureS1A_threshold_stage_profiles.tsv','A'); calls=load('Supplementary_FigureS1B_top_stage_calls.tsv','B'); cont=load('Supplementary_FigureS1C_threshold_contrasts.tsv','C'); counts=load('Supplementary_FigureS1D_gene_counts.tsv','D')
    title(f,'A','Threshold-dependent profiles',.055,.95)
    for j,p in enumerate(MODS):
        a=ax(f,[.105+j*.294,.73,.245,.18]); profile(a,d[d.program==p],thresholds=(50,100,200,500),labels=j==0); a.set_title(module_label(p),color=color(p),weight='bold')
    f.legend([Line2D([],[],color=INK,ls=ls,marker=mk,mfc=INK if mk=='s' else 'white') for ls,mk in [(':','v'),('-.','^'),('--','o'),('-','s')]],['Top 50','Top 100','Top 200','Top 500'],loc='center',bbox_to_anchor=(.6,.97),ncol=4,frameon=False,fontsize=7)
    title(f,'B','Top-stage calls',.055,.65); title(f,'C','Localization contrasts',.53,.65)
    a=ax(f,[.13,.435,.30,.175]); a.set_xlim(0,4); a.set_ylim(3,0)
    stagecolors=['#9EBCBD','#80A8AE','#567F90','#C3B091','#AF8A66','#8C6249']
    for i,p in enumerate(MODS):
        for j,t in enumerate([50,100,200,500]):
            q=calls[(calls.program==p)&(calls.top_n==t)].iloc[0]; si=STAGES.index(q.top_stage)
            a.add_patch(Rectangle((j,i),1,1,facecolor=stagecolors[si],edgecolor='white',lw=2)); a.text(j+.5,i+.5,SHORT[si],ha='center',va='center',color='white',fontsize=9)
    a.set_xticks(np.arange(4)+.5,['Top 50','Top 100','Top 200','Top 500']); a.set_yticks(np.arange(3)+.5,['NTM1','NTM2','NTM3']); a.tick_params(length=0); a.set_xlabel('Gene selection'); [sp.set_visible(False) for sp in a.spines.values()]
    a=ax(f,[.64,.435,.30,.175]); zero(a); clean(a)
    for p in MODS:
        q=cont[cont.program==p].sort_values('top_n'); a.plot(range(4),q.localization_contrast,'o-',color=color(p),ms=4)
    a.set_xticks(range(4),['Top 50','Top 100','Top 200','Top 500']); a.set_xlabel('Gene selection'); a.set_ylabel('Prenatal − postnatal'); a.set_ylim(-.33,.30)
    title(f,'D','Mapped genes and restart channels',.055,.34)
    a=ax(f,[.13,.08,.79,.215]); a.axis('off'); a.set_xlim(0,12); a.set_ylim(4.2,-.7)
    for j,p in enumerate(MODS):
        x=j*4; a.text(x+.1,-.5,module_label(p),color=color(p),weight='bold',fontsize=8)
        for xx,lab in zip([.1,.9,1.9,2.85],['Top','Mapped','Increased','Decreased']): a.text(x+xx,0,lab,fontsize=7,weight='bold')
        for i,t in enumerate([50,100,200,500]):
            g=counts[(counts.program==p)&(counts.top_n==t)]; vals=[t,int(g.n_mapped_to_gene_network.iloc[0]),int(g[g.channel=='positive'].n_restart_genes.iloc[0]),int(g[g.channel=='negative'].n_restart_genes.iloc[0])]
            for xx,v in zip([.1,.9,1.9,2.85],vals): a.text(x+xx,i*.7+.85,str(v),fontsize=8)
            a.plot([x,x+3.7],[i*.7+1.05]*2,color=LIGHT,lw=.6)
    note(f,STAGE_KEY,.08,.035,6.3); export(f,name)

def s4():
    name='Supplementary_Figure_S4'; f=fig(name,8.8)
    pos=load('Supplementary_FigureS4A_q_pos_profiles.tsv','A'); neg=load('Supplementary_FigureS4B_q_neg_profiles.tsv','B'); con=load('Supplementary_FigureS4C_channel_contrasts.tsv','C'); comp=load('Supplementary_FigureS4D_channel_composition.tsv','D')
    for i,(d,label) in enumerate([(pos,'Increased channel (q_pos)'),(neg,'Decreased channel (q_neg)')]):
        y=.745-i*.27; title(f,'AB'[i],label,.055,y+.20)
        for j,p in enumerate(MODS):
            a=ax(f,[.105+j*.294,y,.245,.165]); profile(a,d[d.program==p],'channel_cosine',labels=j==0); a.set_title(module_label(p),weight='bold',color=color(p))
    title(f,'C','Channel contrasts',.055,.39); title(f,'D','Channel composition',.55,.39)
    a=ax(f,[.18,.12,.30,.225]); zero(a,True); clean(a,'x')
    rows=[(p,t) for p in MODS for t in [200,500]]
    for i,(p,t) in enumerate(rows):
        for j,ch in enumerate(['positive','negative']):
            q=con[(con.program==p)&(con.top_n==t)&(con.channel==ch)].iloc[0]
            if pd.notna(q.localization_contrast): a.plot(q.localization_contrast,i+(j-.5)*.20,'o',color=CC[ch],ms=4)
            else: a.text(-.32,i+(j-.5)*.20,'—',fontsize=7,color=GRAY,va='center')
    a.set_yticks(range(6),[f'{mod(p)} · Top {t}' for p,t in rows]); a.set_ylim(5.6,-.6); a.set_xlim(-.34,.34); a.set_xlabel('Prenatal − postnatal contrast')
    a=ax(f,[.68,.12,.27,.225]); channel_composition(a,comp)
    threshold_legend(f,.965)
    f.legend([Line2D([],[],color=CC[x],marker='o') for x in CC],['Increased channel (q_pos)','Decreased channel (q_neg)'],loc='center',bbox_to_anchor=(.5,.055),ncol=2,frameon=False)
    note(f,STAGE_KEY,.08,.019,6.3); export(f,name)

def s6():
    name='Supplementary_Figure_S6'; f=fig(name,7.9)
    data=[load('Supplementary_FigureS6A_marker_PC_adjustment.tsv','A'),load('Supplementary_FigureS6B_all_marker_adjustment.tsv','B')]
    nn=load('Supplementary_FigureS6C_NNLS_adjustment.tsv','C'); cells=load('Supplementary_FigureS6D_celltype_localization.tsv','D')
    limits=(min(float(d.log2_min.min()) for d in data)-.15,max(float(d.log2_max.max()) for d in data)+.15)
    for j,(d,label) in enumerate(zip(data,['Marker-PC adjustment','All-marker adjustment'])):
        title(f,'AB'[j],label,.055+j*.31,.94); a=ax(f,[.11+j*.31,.695,.23,.19])
        for i,q in enumerate(d.itertuples()):
            a.plot([q.log2_min,q.log2_max],[i,i],color=color(q.module),lw=2); a.plot(q.log2_ratio,i,'o',color=color(q.module),ms=5)
        a.set_yticks(range(3),['NTM1','NTM2','NTM3'] if j==0 else ['']*3); a.set_ylim(2.5,-.5); a.set_xlim(*limits); a.set_xticks([-3,-2,-1,0,1]); zero(a,True); clean(a,'x'); a.set_xlabel(r'$\log_2$(adjusted / original β)')
    title(f,'C','NNLS adjustment',.675,.94); a=ax(f,[.745,.695,.21,.19])
    for i,p in enumerate(MODS):
        for j,ad in enumerate(nn.adjustment.unique()):
            q=nn[(nn.program==p)&(nn.adjustment==ad)].iloc[0]; a.plot(q.beta_ratio,i+(j-.5)*.23,'o' if j==0 else 's',color=color(p),mfc='white' if j==0 else color(p),ms=4.5)
    a.set_yticks(range(3),['NTM1','NTM2','NTM3']); a.set_ylim(2.5,-.5); a.set_xlim(0,1.1); a.axvline(1,color=GRAY,ls='--',lw=.7); clean(a,'x'); a.set_xlabel('Adjusted / original β')
    note(f,'A–B: source estimate and cohort min–max range (not CI).   C: ○ NNLS-PC   ■ NNLS-proportion.',.08,.63,6.8)
    title(f,'D','Cell-type context',.055,.575)
    for j,ds in enumerate(cells.dataset.unique()):
        g=cells[cells.dataset==ds]; ct=list(g.celltype.unique()); mods=list(g.module.unique()); cols=[(m,idx) for m in mods for idx in sorted(g.source_profile_index.unique())]
        profile_labels={1:'Top 200',2:'Top 500'}
        vals=[[g[(g.celltype==c)&(g.module==m)&(g.source_profile_index==idx)].value.iloc[0] for m,idx in cols] for c in ct]
        a=ax(f,[.22+j*.45,.235,.28,.28]); im=heat(a,vals,[c.replace('Excitatory neurons','Excitatory').replace('Inhibitory neurons','Inhibitory') for c in ct],[f'{mod(m)}\n{profile_labels[idx]}' for m,idx in cols],-1.6,1.6,DIV,True,6.3)
        a.tick_params(axis='y',labelsize=7); a.tick_params(axis='x',labelsize=7); a.set_title(ds,loc='left',pad=10,weight='bold')
        for v in [2,4]: a.axvline(v,color='white',lw=3)
        donor_note = 'Donors: 60 across displayed cell types' if ds == 'PsychENCODE' else 'Donors: microglia 29; other displayed cell types 31'
        note(f,donor_note,.22+j*.45,.185,7)
    cb(f,im,[.34,.115,.38,.011],'Source cell-type localization value')
    note(f,'NTM1 upregulated · NTM2 downregulated · NTM3 mixed',.08,.020,7)
    note(f,'Native NTM-score context analyses. Columns are native NTM gene-set profiles: Top 200 and Top 500.',.08,.045,7)
    export(f,name)
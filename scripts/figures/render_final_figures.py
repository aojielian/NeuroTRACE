#!/usr/bin/env python3
"""Deterministic final figure renderer for supplied panel data."""
from pathlib import Path
import os, hashlib
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/neurotrace_public_mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from PIL import Image, ImageDraw

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'data'/'processed'/'plotdata'
OUT=ROOT/'figures'
MODS=['NTM1_ASD_up','NTM2_ASD_down','NTM3_ASD_signed']
MC=['#287A80','#B65E36','#7861A4']
CC={'positive':'#B55543','negative':'#336EAA'}
STAGES=['early_prenatal','mid_prenatal','late_prenatal','childhood','adolescence','adulthood']
SL=['Early pre.','Mid pre.','Late pre.','Childhood','Adolesc.','Adult']
SHORT=['EP','MP','LP','CH','AD','AU']
STAGE_KEY='EP early prenatal   MP mid prenatal   LP late prenatal   CH childhood   AD adolescence   AU adulthood'
INK='#27313A'; GRAY='#7D858C'; LIGHT='#E8ECEE'
CMAP=LinearSegmentedColormap.from_list('teal', ['#F4F7F6','#A3C8C5','#17636A'])
DIV=LinearSegmentedColormap.from_list('signed',['#356A9B','#F8F8F6','#B85D42'])
CMAP.set_bad('#E4E6E8'); DIV.set_bad('#E4E6E8')
METHODS=['DIRECT_NATIVE_PROJECTION','DUAL_HEAD_GENE_FIRST_PPR','GENE_ONLY_PPR_BASELINE','GRAPH_OT','HETERO_PPR_STANDARD','SIMULATION_PROXY','MEAN_SIGNATURE','RANDOM','SIGNED_GENE_FIRST_PPR_SINGLE_HEAD','SIMPLE_OT','STAGE_CORRELATION','UNSIGNED_GENE_PPR']
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
    c.outline.set_visible(False); c.ax.tick_params(length=2,labelsize=6.5,pad=2); c.set_label(label,size=7,labelpad=3)
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

def export(f,name):
    folder=OUT/('supplementary' if name.startswith('Supplementary') else 'main')
    folder.mkdir(parents=True,exist_ok=True)
    f.canvas.draw()
    for t in f.findobj(matplotlib.text.Text):
        if t.get_text(): TEXTS.append({'figure':name,'text':t.get_text()})
    pdf=folder/(name+'.pdf'); png=folder/(name+'.png'); tif=folder/(name+'.tiff')
    f.savefig(pdf,metadata={'Creator':'Scientific plotting implementation','Title':name})
    f.savefig(png,dpi=300)
    with Image.open(png) as im:
        im.convert('RGB').save(tif,dpi=(300,300),compression='tiff_lzw')
        w,h=im.size
    for path in [pdf,png,tif]: EXPORTS.append(dict(figure=name,path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),width_px=w if path.suffix!='.pdf' else '',height_px=h if path.suffix!='.pdf' else '',dpi=300 if path.suffix!='.pdf' else 'vector'))
    plt.close(f); print('RENDERED',name,flush=True)

def figure1():
    f=fig('Figure1',5.8)
    d=load('Figure1_method_components.tsv','A–C')
    for suffix in 'ABC': load(f'Figure1{suffix}_method_components.tsv',suffix)
    defs={r.component:r.definition for r in d.itertuples()}
    a=ax(f,[.04,.06,.92,.89]); a.set_xlim(0,100); a.set_ylim(0,100); a.axis('off')
    def tx(x,y,s,sz=9,**kw): return a.text(x,y,s,fontsize=sz,ha='center',va='center',**kw)
    def arrow(x,y,x2,y2,c=GRAY): a.add_patch(FancyArrowPatch((x,y),(x2,y2),arrowstyle='-|>',mutation_scale=10,lw=1.1,color=c))
    tx(0,98,'A',11,weight='bold',ha='left') if False else None
    title(f,'A','Signed restarts',.045,.94); title(f,'B','Gene-only diffusion',.49,.94)
    tx(14,83,'Signed native module weights',9,weight='bold'); tx(14,75,r'$w_g$',15)
    arrow(14,70,5,60,CC['positive']); arrow(15,70,30,60,CC['negative'])
    tx(7,56,'Positive restart',9,color=CC['positive']); tx(30,56,'Negative restart',9,color=CC['negative'])
    tx(7,49,r'$r^+_g\propto |w_g|,\ w_g>0$',10); tx(30,49,r'$r^-_g\propto |w_g|,\ w_g<0$',10)
    tx(18,39,'Normalize each non-empty channel',8)
    # Diagram-only gene topology; no empirical gene identities or numerical values.
    pts=np.array([[56,75],[63,86],[72,83],[79,72],[70,64],[60,64],[67,74],[84,83]])
    edges=[(0,1),(0,5),(0,6),(1,2),(1,6),(2,3),(2,6),(2,7),(3,4),(3,7),(4,5),(4,6),(5,6)]
    for i,j in edges: a.plot(pts[[i,j],0],pts[[i,j],1],color='#BCC5CB',lw=1,zorder=1)
    for x,y in pts: a.add_patch(Circle((x,y),1.55,facecolor='white',edgecolor=INK,lw=.9,zorder=2))
    tx(71,94,'Gene similarity network',9); tx(86,63,'Schematic',6.5,color=GRAY)
    arrow(43,58,55,73,CC['positive']); arrow(43,48,57,64,CC['negative'])
    arrow(67,61,64,49,CC['positive']); arrow(73,61,82,49,CC['negative'])
    tx(63,44,r'$q_{\mathrm{pos}}$',13,color=CC['positive']); tx(84,44,r'$q_{\mathrm{neg}}$',13,color=CC['negative'])
    a.plot([0,100],[35,35],color=LIGHT,lw=1)
    title(f,'C','Dual-head readouts',.045,.335)
    tx(21,24,r'$q_{\mathrm{magnitude}} = q_{\mathrm{pos}} + q_{\mathrm{neg}}$',12,weight='bold')
    tx(72,24,r'$q_{\mathrm{signed}} = q_{\mathrm{pos}} - q_{\mathrm{neg}}$',12,weight='bold')
    arrow(21,20,21,14); arrow(72,20,60,14); arrow(73,20,88,14)
    tx(21,10,'Developmental localization',9,weight='bold')
    tx(21,3,r'$\cos(q_{\mathrm{magnitude}},p_t)$'+'  ·  six stage profiles',9)
    tx(59,10,'Gene prioritization',9,weight='bold'); tx(59,3,r'$|q_{\mathrm{signed}}(g)|$',10)
    tx(88,10,'Sample scoring',9,weight='bold'); tx(88,3,'Signed gene weights',8.5)
    export(f,'Figure1')

def figure2():
    f=fig('Figure2',10.3)
    d=load('Figure2A_condition_matrix.tsv','A'); dev=load('Figure2B_developmental_metrics.tsv','B'); gene=load('Figure2C_gene_metrics.tsv','C'); tr=load('Figure2D_tradeoff.tsv','D')
    title(f,'A','Simulation conditions',.05,.955)
    dims={'batch_sd':'Batch SD','cohort_noise_sd':'Cohort noise','composition_sd':'Comp. SD','dropout_rate':'Dropout','false_prior_rate':'False prior','adversarial_composition':'Adv. comp.','signal_strength':'Signal','state_overlap':'Overlap','true_prior_recall':'Prior recall'}
    for k,(b,x,w) in enumerate([('step02A',.18,.29),('step02B',.65,.30)]):
        g=d[d.benchmark==b]; cols=list(dict.fromkeys(g.perturbation_dimension)); rows=list(dict.fromkeys(g.condition))
        raw=g.pivot(index='condition',columns='perturbation_dimension',values='raw_value').loc[rows,cols]
        scaled=g.pivot(index='condition',columns='perturbation_dimension',values='within_design_scaled_value').loc[rows,cols]
        a=ax(f,[x,.755,w,.167]); im=heat(a,scaled,[COND[s] for s in rows],[dims[s] for s in cols],0,1,annot=False)
        for i in range(len(rows)):
            for j in range(len(cols)): a.text(j+.5,i+.5,f'{raw.iloc[i,j]:.2f}',ha='center',va='center',fontsize=6.3,color='white' if scaled.iloc[i,j]>.62 else INK)
        a.tick_params(axis='y',labelsize=6.5); a.tick_params(axis='x',labelsize=6.4)
        plt.setp(a.get_xticklabels(),rotation=45,ha='right')
        a.set_title('High-signal design' if k==0 else 'Challenging design',loc='left',weight='bold',pad=8)
    note(f,'Cells: raw parameters; color: supplied within-design scale (0–1).',.18,.69,7)
    title(f,'B','Developmental performance',.05,.65); title(f,'C','Gene prioritization',.54,.65)
    mat=[]
    for m in METHODS:
        row=[]
        for b in ['step02A','step02B']:
            q=dev[(dev.method==m)&(dev.benchmark==b)].iloc[0]; row += [q.alignment_hit_top1_mean,q.true_stage_mass_mean]
        mat.append(row)
    a=ax(f,[.23,.355,.265,.265]); im=heat(a,mat,ML,['Hit@1','Mass','Hit@1','Mass'],0,1,fs=6.8)
    a.axvline(2,color='white',lw=4); a.text(1,-.6,'High signal',ha='center',fontsize=7); a.text(3,-.6,'Challenging',ha='center',fontsize=7)
    mat=[]
    for m in METHODS:
        row=[]
        for b in ['step02A','step02B']:
            q=gene[(gene.method==m)&(gene.benchmark==b)].iloc[0]; row += [q.gene_auroc_mean,q.gene_auprc_mean,q.precision_at_50_mean]
        mat.append(row)
    a=ax(f,[.565,.355,.385,.265]); im=heat(a,mat,['']*12,['AUROC','AUPRC','P@50']*2,0,1,fs=6.8); a.axvline(3,color='white',lw=4)
    a.text(1.5,-.6,'High signal',ha='center',fontsize=7); a.text(4.5,-.6,'Challenging',ha='center',fontsize=7)
    cb(f,im,[.62,.318,.25,.008],'Performance (0–1)')
    note(f,f'n={int(dev.n.iloc[0]):,} per method and design.\n— metric undefined.',.23,.298,6.7)
    title(f,'D','Challenging-design trade-off',.05,.267)
    a=ax(f,[.10,.065,.42,.17]); clean(a); a.set_xlim(0,.72); a.set_ylim(0,.44)
    a.set_xlabel('Developmental alignment hit@1'); a.set_ylabel('Gene AUPRC')
    offsets={'NP':(-13,12),'MS':(-25,-13),'DH':(5,10),'GP':(5,7),'GO':(5,6),'HP':(8,-1),'SP':(-17,8),'RD':(5,4),'SH':(8,2),'SO':(4,11),'SC':(-17,-12),'UP':(6,-14)}
    for m,code in zip(METHODS,CODES):
        q=tr[tr.method==m].iloc[0]; col=INK if code=='DH' else '#638B90'
        a.scatter(q.alignment_hit_top1_mean,q.gene_auprc_mean,s=34 if code=='DH' else 20,color=col,marker='D' if code=='DH' else 'o',edgecolor='white',lw=.4,zorder=3)
        a.annotate(code,(q.alignment_hit_top1_mean,q.gene_auprc_mean),xytext=offsets[code],textcoords='offset points',fontsize=6.5,arrowprops={'arrowstyle':'-','lw':.4,'color':GRAY})
    for i,(code,label) in enumerate(zip(CODES,ML)):
        f.text(.585+(i//6)*.205,.228-(i%6)*.025,f'{code}  {label}',fontsize=6.6,weight='bold' if code=='DH' else 'normal')
    export(f,'Figure2')

def null_rows(d): return d.sort_values(['program','top_n']).reset_index(drop=True)
def figure3():
    f=fig('Figure3',8.5)
    d=load('Figure3A_stage_profiles.tsv','A'); s=null_rows(load('Figure3B_contrast_null_interval.tsv','B')); r=load('Figure3C_empirical_null_replicates.tsv','C'); freq=null_rows(load('Figure3D_observed_stage_frequency.tsv','D'))
    title(f,'A','Developmental profiles',.055,.95); threshold_legend(f,.965)
    for j,p in enumerate(MODS):
        a=ax(f,[.105+j*.294,.73,.245,.18]); profile(a,d[d.program==p],labels=j==0); a.set_title(mod(p),color=color(p),weight='bold')
    title(f,'B','Contrast and null interval',.055,.645); title(f,'C','Empirical null distributions',.56,.645)
    a=ax(f,[.16,.35,.345,.265]); zero(a,True); clean(a,'x')
    labels=[]
    for i,q in s.iterrows():
        c=color(q.program); a.plot([q.null_localization_contrast_q025,q.null_localization_contrast_q975],[i,i],color='#9AA6AE',lw=5,solid_capstyle='round')
        a.plot(q.null_localization_contrast_median,i,'|',color=INK,ms=9,mew=1)
        a.plot(q.observed_localization_contrast,i,'o' if q.top_n==200 else 's',color=c,mfc='white' if q.top_n==200 else c,ms=5,zorder=3)
        labels.append(f'{mod(q.program)} · {q.top_n}')
    a.set_yticks(range(6),labels); a.set_ylim(5.6,-.6); a.set_xlim(-.34,.40); a.set_xlabel('Prenatal − postnatal contrast')
    a=ax(f,[.64,.35,.305,.265]); clean(a,'x'); a.set_xlim(-.34,.40); a.set_ylim(-.45,6.1)
    bins=np.linspace(-.34,.4,61)
    for i,q in s.iterrows():
        vals=r[(r.program==q.program)&(r.top_n==q.top_n)].localization_contrast
        counts,edges=np.histogram(vals,bins=bins); yy=5-i; heights=counts/max(counts)*.7
        a.stairs(yy+heights,edges,baseline=yy,fill=True,color=color(q.program),alpha=.45,lw=.7)
        a.plot([q.observed_localization_contrast]*2,[yy,yy+.77],color=color(q.program),lw=1.4)
    a.set_yticks(np.arange(6)+.15,labels[::-1]); a.tick_params(axis='y',labelsize=6.8); a.set_xlabel('Prenatal − postnatal contrast')
    note(f,'B: point, observed; thick line, null 95% interval; tick, median.  C: histogram heights scaled within each distribution.',.08,.289,6.6)
    title(f,'D','Observed top-stage frequency',.055,.26)
    a=ax(f,[.17,.09,.46,.14]); clean(a,'x'); a.set_xlim(-.025,1.025); a.set_xticks(np.linspace(0,1,6)); a.set_ylim(5.6,-.6)
    for i,q in freq.iterrows():
        a.plot([0,q.null_top_stage_frequency_observed],[i,i],color=LIGHT,lw=2)
        a.plot(q.null_top_stage_frequency_observed,i,'o',color=color(q.program),ms=4)
    a.set_yticks(range(6),labels); a.set_xlabel('Matched-null frequency of the observed top stage')
    note(f,'Stage · frequency · two-sided contrast FDR',.69,.237,6.3)
    for i,q in freq.iterrows():
        f.text(.69,.211-i*.021,f'{SHORT[STAGES.index(q.observed_top_stage)]}  {q.null_top_stage_frequency_observed:.3f}   FDR {fdr(q.BH_FDR_two_sided)}',fontsize=6.9)
    note(f,STAGE_KEY,.08,.02,6.3)
    export(f,'Figure3')

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
    a.set_yticks(range(6),[f'{mod(p)} · {t}' for p,t in rows]); a.set_ylim(5.6,-.6); a.set_xlim(0,1.33)
    a.set_xticks([0,.5,1],['0','0.5','1']); a.set_xlabel('Share of absolute restart weight')
    a.text(1.035,-.72,'Genes + / −',fontsize=6.6); clean(a,'x')

def figure4():
    f=fig('Figure4',8.4)
    d=load('Figure4A_NTM3_channel_profiles.tsv','A'); h=load('Figure4B_all_channel_heatmap.tsv','B'); c=load('Figure4C_readout_contrasts.tsv','C'); comp=load('Figure4D_channel_composition.tsv','D')
    title(f,'A','NTM3 channel profiles',.055,.95); title(f,'B','All-module channel structure',.55,.95)
    a=ax(f,[.10,.665,.355,.245]); profile(a,d,'channel_cosine','channel')
    a.legend([Line2D([],[],color=CC[x]) for x in CC],['Positive','Negative'],loc='upper right',frameon=False,fontsize=7)
    rows=[(p,t,ch) for p in MODS for t in [200,500] for ch in ['positive','negative']]
    mat=[h[(h.program==p)&(h.top_n==t)&(h.channel==ch)].set_index('state_id').loc[STAGES,'channel_cosine'].values for p,t,ch in rows]
    a=ax(f,[.665,.665,.285,.245]); im=heat(a,mat,[f'{mod(p)} {t} {"+" if ch=="positive" else "−"}' for p,t,ch in rows],SHORT,-.32,.32,DIV,False)
    a.tick_params(axis='y',labelsize=6.2)
    for i,row in enumerate(mat):
        if np.isnan(row).all(): a.text(3,i+.5,'Absent',ha='center',va='center',fontsize=6,color=GRAY)
    cb(f,im,[.69,.61,.23,.009],'Channel cosine')
    title(f,'C','Readout contrasts',.055,.545); title(f,'D','Channel composition',.55,.545)
    a=ax(f,[.15,.11,.335,.38]); zero(a,True); clean(a,'x')
    readouts=['q_pos','q_neg','q_magnitude','q_signed']; rl=['Positive','Negative','Magnitude','Signed orientation']
    for j,ro in enumerate(readouts):
        for k,p in enumerate(MODS):
            for t,off,mk in [(200,-.105,'o'),(500,.105,'s')]:
                q=c[(c.program==p)&(c.top_n==t)&(c.readout==ro)].iloc[0]; y=j*3.9+k+off
                if pd.notna(q.localization_contrast): a.plot(q.localization_contrast,y,mk,color=color(p),mfc='white' if t==200 else color(p),ms=4)
                else: a.text(-.325,y,'—',fontsize=7,color=GRAY,va='center')
        if j<3: a.axhline(j*3.9+2.95,color=LIGHT,lw=.7)
    a.set_yticks([j*3.9+1 for j in range(4)],rl); a.tick_params(axis='y',labelsize=7); a.set_ylim(14.6,-.9); a.set_xlim(-.35,.35); a.set_xlabel('Prenatal − postnatal contrast')
    a=ax(f,[.67,.20,.28,.285]); channel_composition(a,comp)
    note(f,'q_signed: signed directional orientation',.08,.051,7)
    threshold_legend(f,.57); module_legend(f,.038)
    note(f,STAGE_KEY,.08,.012,6.3)
    export(f,'Figure4')

def forest(a,d,cohort=None,allrows=False,limits=None):
    if cohort: d=d[d.cohort==cohort]
    if allrows:
        groups=[(co,p,t) for co in d.cohort.unique() for p in MODS for t in [200,500]]
    else: groups=[(cohort,p,t) for p in MODS for t in [200,500]]
    labels=[]
    for i,(co,p,t) in enumerate(groups):
        g=d[(d.cohort==co)&(d.program==p)&(d.top_n==t)]
        for j,q in enumerate(g.itertuples()):
            y=i+(j-(len(g)-1)/2)*.25
            a.plot([q.ci_low,q.ci_high],[y,y],color=color(p),lw=1.1)
            a.plot(q.standardized_beta,y,'o' if j==0 else 's',color=color(p),mfc='white' if j==0 and len(g)>1 else color(p),ms=4)
        labels.append(f'{mod(p)} · {t}')
    if limits is None: limits=(min(-.9,float(d.ci_low.min())-.12),max(2.35,float(d.ci_high.max())+.12))
    a.set_yticks(range(len(groups)),labels); a.set_ylim(len(groups)-.45,-.65); a.set_xlim(*limits); a.set_xlabel('Standardized β (95% CI)'); zero(a,True); clean(a,'x')
    if allrows: a.axhline(5.5,color=GRAY,lw=.6)

def figure5():
    f=fig('Figure5',8.6)
    d=load('Figure5A_forest.tsv','A'); pair=load('Figure5B_paired_beta.tsv','B'); matrix=load('Figure5C_effect_FDR_matrix.tsv','C'); su=load('Figure5D_cohort_summary.tsv','D')
    title(f,'A','External standardized effects',.055,.95)
    for j,co in enumerate(d.cohort.unique()):
        a=ax(f,[.13+j*.46,.68,.34,.225]); forest(a,d,co)
        q=d[d.cohort==co].iloc[0]; a.set_title(f'{co}   n={q.n} ({q.n_ASD} ASD / {q.n_Control} control)',fontsize=8,loc='left',pad=9)
    f.legend([Line2D([],[],color=INK,marker='o',mfc='white'),Line2D([],[],color=INK,marker='s')],['Native signed module','Signed gene PPR'],loc='center',bbox_to_anchor=(.52,.623),ncol=2,frameon=False)
    title(f,'B','Native versus signed PPR',.055,.572); title(f,'C','Effects and FDR',.54,.572)
    a=ax(f,[.11,.315,.32,.215]); a.plot([0,1.6],[0,1.6],color=GRAY,ls='--',lw=.8)
    for q in pair.itertuples():
        a.plot(q.NATIVE_SIGNED_MODULE_standardized_beta,q.SIGNED_GENE_PPR_standardized_beta,'o' if q.top_n==200 else 's',color=color(q.program),mfc='white' if q.top_n==200 else color(q.program),ms=4.5)
    a.set_xlim(0,1.6); a.set_ylim(0,1.6); a.set_aspect('equal'); a.set_xlabel('Native standardized β'); a.set_ylabel('Signed-PPR standardized β'); clean(a)
    note(f,'B: ○ top 200   ■ top 500',.125,.263,6.7)
    rows=[(p,t) for p in MODS for t in [200,500]]; cols=[(co,m) for co in matrix.cohort.unique() for m in matrix.method.unique()]
    vals=[]; qs=[]
    for p,t in rows:
        gr=[]; fq=[]
        for co,m in cols:
            q=matrix[(matrix.program==p)&(matrix.top_n==t)&(matrix.cohort==co)&(matrix.method==m)].iloc[0]; gr.append(q.standardized_beta); fq.append(q.fdr)
        vals.append(gr); qs.append(fq)
    a=ax(f,[.65,.315,.29,.215]); im=heat(a,vals,[f'{mod(p)} {t}' for p,t in rows],['Native','PPR']*2,0,1.6,annot=False)
    for i in range(6):
        for j in range(4): a.text(j+.5,i+.5,f'{vals[i][j]:.2f}'+('*' if qs[i][j]<.05 else ''),ha='center',va='center',fontsize=7,color='white' if vals[i][j]>.95 else INK)
    a.axvline(2,color='white',lw=4); a.text(1,-.5,'GSE102741',ha='center',fontsize=7); a.text(3,-.5,'GSE64018',ha='center',fontsize=7)
    note(f,'Cell: standardized β; * BH-FDR < 0.05 (36-model family).',.59,.276,6.5)
    title(f,'D','Cohort-level summary',.055,.232)
    a=ax(f,[.23,.07,.66,.125]); fields=['direction_concordant_n','native_fdr_significant_n','signed_ppr_fdr_significant_n']; labs=['Direction concordant','Native FDR < 0.05','PPR FDR < 0.05']
    for j,q in enumerate(su.itertuples()):
        for i,key in enumerate(fields):
            v=getattr(q,key); y=i+(j-.5)*.23; a.plot(v,y,'o' if j==0 else 's',color=['#537E8B','#BA7044'][j],ms=5); a.text(v+.13,y,f'{v}/{q.comparison_denominator}',va='center',fontsize=7)
    a.set_yticks(range(3),labs); a.set_ylim(2.6,-.6); a.set_xlim(-.1,7); a.set_xticks(range(7)); a.set_xlabel('Module–threshold comparisons'); clean(a,'x')
    f.legend([Line2D([],[],color='#537E8B',marker='o',ls=''),Line2D([],[],color='#BA7044',marker='s',ls='')],list(su.cohort),loc='center',bbox_to_anchor=(.63,.235),ncol=2,frameon=False,fontsize=7)
    module_legend(f,.018); export(f,'Figure5')

def s1():
    name='Supplementary_Figure_S1'; f=fig(name,8.3)
    d=load('Supplementary_FigureS1A_threshold_stage_profiles.tsv','A'); calls=load('Supplementary_FigureS1B_top_stage_calls.tsv','B'); cont=load('Supplementary_FigureS1C_threshold_contrasts.tsv','C'); counts=load('Supplementary_FigureS1D_gene_counts.tsv','D')
    title(f,'A','Threshold-dependent profiles',.055,.95)
    for j,p in enumerate(MODS):
        a=ax(f,[.105+j*.294,.73,.245,.18]); profile(a,d[d.program==p],thresholds=(50,100,200,500),labels=j==0); a.set_title(mod(p),color=color(p),weight='bold')
    f.legend([Line2D([],[],color=INK,ls=ls,marker=mk,mfc=INK if mk=='s' else 'white') for ls,mk in [(':','v'),('-.','^'),('--','o'),('-','s')]],['Top 50','Top 100','Top 200','Top 500'],loc='center',bbox_to_anchor=(.6,.97),ncol=4,frameon=False,fontsize=7)
    title(f,'B','Top-stage calls',.055,.65); title(f,'C','Localization contrasts',.53,.65)
    a=ax(f,[.13,.435,.30,.175]); a.set_xlim(0,4); a.set_ylim(3,0)
    stagecolors=['#9EBCBD','#80A8AE','#567F90','#C3B091','#AF8A66','#8C6249']
    for i,p in enumerate(MODS):
        for j,t in enumerate([50,100,200,500]):
            q=calls[(calls.program==p)&(calls.top_n==t)].iloc[0]; si=STAGES.index(q.top_stage)
            a.add_patch(Rectangle((j,i),1,1,facecolor=stagecolors[si],edgecolor='white',lw=2)); a.text(j+.5,i+.5,SHORT[si],ha='center',va='center',color='white',fontsize=9)
    a.set_xticks(np.arange(4)+.5,['50','100','200','500']); a.set_yticks(np.arange(3)+.5,['NTM1','NTM2','NTM3']); a.tick_params(length=0); a.set_xlabel('Top genes'); [sp.set_visible(False) for sp in a.spines.values()]
    a=ax(f,[.64,.435,.30,.175]); zero(a); clean(a)
    for p in MODS:
        q=cont[cont.program==p].sort_values('top_n'); a.plot(range(4),q.localization_contrast,'o-',color=color(p),ms=4)
    a.set_xticks(range(4),['50','100','200','500']); a.set_xlabel('Top genes'); a.set_ylabel('Prenatal − postnatal'); a.set_ylim(-.33,.30)
    title(f,'D','Mapped genes and restart channels',.055,.34)
    a=ax(f,[.13,.08,.79,.215]); a.axis('off'); a.set_xlim(0,12); a.set_ylim(4.2,-.7)
    for j,p in enumerate(MODS):
        x=j*4; a.text(x+.1,-.5,mod(p),color=color(p),weight='bold',fontsize=8)
        for xx,lab in zip([.1,.9,1.9,2.85],['Top','Mapped','Positive','Negative']): a.text(x+xx,0,lab,fontsize=7,weight='bold')
        for i,t in enumerate([50,100,200,500]):
            g=counts[(counts.program==p)&(counts.top_n==t)]; vals=[t,int(g.n_mapped_to_gene_network.iloc[0]),int(g[g.channel=='positive'].n_restart_genes.iloc[0]),int(g[g.channel=='negative'].n_restart_genes.iloc[0])]
            for xx,v in zip([.1,.9,1.9,2.85],vals): a.text(x+xx,i*.7+.85,str(v),fontsize=8)
            a.plot([x,x+3.7],[i*.7+1.05]*2,color=LIGHT,lw=.6)
    note(f,STAGE_KEY,.08,.035,6.3); export(f,name)

def s2():
    name='Supplementary_Figure_S2'; f=fig(name,9.5)
    specs=[('A','alignment_hit_top1_mean','Developmental alignment'),('B','true_stage_mass_mean','True-stage mass'),('C','gene_auprc_mean','Gene AUPRC'),('D','precision_at_50_mean','Precision at 50')]
    for k,(letter,metric,label) in enumerate(specs):
        d=load(f'Supplementary_FigureS2{letter}_{metric}.tsv',letter)
        rows=[]
        for b in ['step02A','step02B']:
            for cond in d[d.benchmark==b].condition.unique(): rows.append((b,cond))
        mat=[[d[(d.benchmark==b)&(d.condition==c)&(d.method==m)][metric].iloc[0] for m in METHODS] for b,c in rows]
        col=k%2; row=k//2; x=.19+col*.455; y=.625-row*.40
        title(f,letter,label,.055+col*.455,y+.30)
        a=ax(f,[x,y,.315,.27]); im=heat(a,mat,[('H · ' if b=='step02A' else 'C · ')+COND[c] for b,c in rows],CODES,0,1,annot=False)
        a.tick_params(axis='y',labelsize=6.3); a.tick_params(axis='x',labelsize=6.2,pad=5); a.axhline(6,color='white',lw=4)
        for i in range(12):
            for j in range(12):
                if pd.isna(mat[i][j]): a.text(j+.5,i+.5,'—',ha='center',va='center',fontsize=6,color=GRAY)
    cb(f,im,[.34,.135,.39,.01],'Performance (0–1); gray / — = undefined')
    for i,(code,label) in enumerate(zip(CODES,ML)): f.text(.08+(i//4)*.31,.095-(i%4)*.02,f'{code}  {label}',fontsize=6.8)
    note(f,f'H high-signal design   C challenging design   Each cell: supplied mean of {int(d.n_replicates.iloc[0])} replicates.',.08,.015,6.8)
    export(f,name)

def s3():
    name='Supplementary_Figure_S3'; f=fig(name,8.6)
    r=load('Supplementary_FigureS3_empirical_null_distributions.tsv','A–F'); s=load('Figure3B_contrast_null_interval.tsv','A–F summary')
    bins=np.linspace(-.34,.40,61)
    for i,p in enumerate(MODS):
        for j,t in enumerate([200,500]):
            k=i*2+j; x=.11+j*.465; y=.725-i*.285
            title(f,chr(65+k),f'{mod(p)} · top {t}',.055+j*.465,y+.21)
            a=ax(f,[x,y,.365,.18]); q=s[(s.program==p)&(s.top_n==t)].iloc[0]; vals=r[(r.program==p)&(r.top_n==t)].localization_contrast
            a.hist(vals,bins=bins,color=color(p),alpha=.65,edgecolor='white',lw=.25)
            a.axvspan(q.null_localization_contrast_q025,q.null_localization_contrast_q975,color=color(p),alpha=.10)
            a.axvline(q.null_localization_contrast_q025,color=GRAY,ls=':',lw=.85); a.axvline(q.null_localization_contrast_q975,color=GRAY,ls=':',lw=.85)
            a.axvline(q.null_localization_contrast_median,color=INK,ls='--',lw=.9); a.axvline(q.observed_localization_contrast,color=color(p),lw=1.7)
            a.set_xlim(-.34,.4); clean(a); a.set_ylabel('Null replicates'); a.set_xlabel('Prenatal − postnatal contrast')
            a.text(.025,1.025,f'n={q.n_null_replicates:,}   Two-sided FDR={fdr(q.BH_FDR_two_sided)}',transform=a.transAxes,va='bottom',fontsize=6.5)
    f.legend([Line2D([],[],color=INK,lw=1.7),Line2D([],[],color=INK,ls='--'),Line2D([],[],color=GRAY,ls=':')],['Observed','Null median','Null 2.5th / 97.5th percentiles'],loc='center',bbox_to_anchor=(.52,.045),ncol=3,frameon=False)
    export(f,name)

def s4():
    name='Supplementary_Figure_S4'; f=fig(name,8.8)
    pos=load('Supplementary_FigureS4A_q_pos_profiles.tsv','A'); neg=load('Supplementary_FigureS4B_q_neg_profiles.tsv','B'); con=load('Supplementary_FigureS4C_channel_contrasts.tsv','C'); comp=load('Supplementary_FigureS4D_channel_composition.tsv','D')
    for i,(d,label) in enumerate([(pos,'Positive-channel profiles'),(neg,'Negative-channel profiles')]):
        y=.745-i*.27; title(f,'AB'[i],label,.055,y+.20)
        for j,p in enumerate(MODS):
            a=ax(f,[.105+j*.294,y,.245,.165]); profile(a,d[d.program==p],'channel_cosine',labels=j==0); a.set_title(mod(p),weight='bold',color=color(p))
    title(f,'C','Channel contrasts',.055,.39); title(f,'D','Channel composition',.55,.39)
    a=ax(f,[.18,.12,.30,.225]); zero(a,True); clean(a,'x')
    rows=[(p,t) for p in MODS for t in [200,500]]
    for i,(p,t) in enumerate(rows):
        for j,ch in enumerate(['positive','negative']):
            q=con[(con.program==p)&(con.top_n==t)&(con.channel==ch)].iloc[0]
            if pd.notna(q.localization_contrast): a.plot(q.localization_contrast,i+(j-.5)*.20,'o',color=CC[ch],ms=4)
            else: a.text(-.32,i+(j-.5)*.20,'—',fontsize=7,color=GRAY,va='center')
    a.set_yticks(range(6),[f'{mod(p)} · {t}' for p,t in rows]); a.set_ylim(5.6,-.6); a.set_xlim(-.34,.34); a.set_xlabel('Prenatal − postnatal contrast')
    a=ax(f,[.68,.12,.27,.225]); channel_composition(a,comp)
    threshold_legend(f,.965)
    f.legend([Line2D([],[],color=CC[x],marker='o') for x in CC],['Positive channel','Negative channel'],loc='center',bbox_to_anchor=(.5,.055),ncol=2,frameon=False)
    note(f,STAGE_KEY,.08,.019,6.3); export(f,name)

def s5():
    name='Supplementary_Figure_S5'; f=fig(name,9.1)
    methods=['NATIVE_SIGNED_MODULE','SIGNED_GENE_PPR','UNSIGNED_GENE_PPR']; names=['Native signed module','Signed gene PPR','Unsigned gene PPR']
    all_d=[load(f'Supplementary_FigureS5{"ABC"[j]}_{m}.tsv','ABC'[j]) for j,m in enumerate(methods)]
    limits=(min(float(d.ci_low.min()) for d in all_d)-.12,max(float(d.ci_high.max()) for d in all_d)+.12)
    for j,m in enumerate(methods):
        d=all_d[j]
        x=.11+j*.295; title(f,'ABC'[j],names[j],.055+j*.295,.95)
        a=ax(f,[x,.55,.245,.355]); forest(a,d,allrows=True,limits=limits); a.tick_params(axis='y',labelsize=6.5); a.tick_params(axis='x',labelsize=6.5)
        if j>0: a.set_yticklabels([])
        a.text(.02,1.035,'GSE102741',transform=a.transAxes,va='bottom',fontsize=7,weight='bold'); a.text(.02,.496,'GSE64018',transform=a.transAxes,va='center',fontsize=7,weight='bold',bbox={'facecolor':'white','edgecolor':'none','pad':1})
    title(f,'D','Gene mapping and score dispersion',.055,.465)
    d=load('Supplementary_FigureS5D_mapping_and_score_sd.tsv','D')
    rows=[(co,p,t) for co in d.cohort.unique() for p in MODS for t in [200,500]]
    for k,(val,label,x) in enumerate([('n_common_genes','Mapped genes',.22),('score_SD_before_standardization','Score SD before standardization',.65)]):
        mat=[[d[(d.cohort==co)&(d.program==p)&(d.top_n==t)&(d.method==m)][val].iloc[0] for m in methods] for co,p,t in rows]
        a=ax(f,[x,.095,.29,.315]); im=heat(a,mat,[f'{"102741" if co=="GSE102741" else "64018"} · {mod(p)} {t}' for co,p,t in rows] if k==0 else ['']*12,['Native','Signed PPR','Unsigned PPR'],0,None,annot=False)
        for i in range(12):
            for j in range(3): a.text(j+.5,i+.5,f'{int(mat[i][j]):,}' if k==0 else f'{mat[i][j]:.3f}',ha='center',va='center',fontsize=7,color='white' if im.norm(mat[i][j])>.62 else INK)
        a.axhline(6,color='white',lw=4); a.set_title(label,loc='left',pad=8,fontsize=8); a.tick_params(axis='y',labelsize=6.5)
    note(f,'All supplied models and both thresholds; intervals are source 95% CIs.',.08,.035,7)
    module_legend(f,.016); export(f,name)

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
        vals=[[g[(g.celltype==c)&(g.module==m)&(g.source_profile_index==idx)].value.iloc[0] for m,idx in cols] for c in ct]
        a=ax(f,[.22+j*.45,.235,.28,.28]); im=heat(a,vals,[c.replace('Excitatory neurons','Excitatory').replace('Inhibitory neurons','Inhibitory') for c in ct],[f'{mod(m)}·{idx}' for m,idx in cols],-1.6,1.6,DIV,True,6.3)
        a.tick_params(axis='y',labelsize=7); a.tick_params(axis='x',labelsize=6); a.set_title(ds,loc='left',pad=10,weight='bold')
        for v in [2,4]: a.axvline(v,color='white',lw=3)
        donors=sorted(g.n_donors.unique()); note(f,'Donors: '+', '.join(str(x) for x in donors),.22+j*.45,.185,7)
    cb(f,im,[.34,.115,.38,.011],'Source cell-type localization value')
    note(f,'Native NTM-score context analyses. Columns 1/2 retain the two source profiles; indices are not biological groups.',.08,.045,6.7)
    export(f,name)

def s7():
    name='Supplementary_Figure_S7'; f=fig(name,6.7)
    d=load('Supplementary_FigureS7A_cross_disorder_heatmap.tsv','A'); su=load('Supplementary_FigureS7B_disease_summary.tsv','B'); md=load('Supplementary_FigureS7C_MDD_forest.tsv','C')
    title(f,'A','Cross-disorder effects',.055,.94); title(f,'B','Disease-level concordance',.53,.94)
    diseases=['SCZ','BD','MDD']; modules=list(d.module.unique()); vals=[[d[(d.disease==di)&(d.module==m)].beta.iloc[0] for m in modules] for di in diseases]
    a=ax(f,[.13,.66,.32,.23]); im=heat(a,vals,diseases,[mod(m) for m in modules],-.35,.35,DIV,False)
    for i,di in enumerate(diseases):
        for j,m in enumerate(modules):
            q=d[(d.disease==di)&(d.module==m)].iloc[0]; lab='' if pd.isna(q.sig_label) else str(q.sig_label)
            a.text(j+.5,i+.5,f'{q.beta:.2f}{lab}',ha='center',va='center',fontsize=9)
    a=ax(f,[.62,.66,.28,.23]); clean(a,'x')
    for i,di in enumerate(diseases):
        q=su[su.disease==di].iloc[0]; a.plot([0,q.concordance],[i,i],color=LIGHT,lw=3); a.plot(q.concordance,i,'o',color='#376F7D',ms=5)
        a.text(q.concordance+.06,i,str(q.label),va='center',fontsize=8)
    a.set_yticks(range(3),diseases); a.set_ylim(2.5,-.5); a.set_xlim(-.04,1.35); a.set_xticks([0,.5,1]); a.set_xlabel('Positive-direction concordance')
    for j,di in enumerate(diseases):
        q=su[su.disease==di].iloc[0]; f.text(.56,.57-j*.028,f'{di}: aggregate FDR {q.stouffer_fdr_across_diseases:.3g}',fontsize=7)
    cb(f,im,[.16,.57,.25,.012],'Native NTM-score β')
    title(f,'C','Available MDD dataset',.055,.435)
    a=ax(f,[.17,.15,.49,.225]); zero(a,True); clean(a,'x')
    for i,q in enumerate(md.itertuples()):
        a.plot([q.ci_low,q.ci_high],[i,i],color=color(q.program),lw=1.5); a.plot(q.beta_target_vs_Control,i,'o',color=color(q.program),ms=5)
        f.text(.73,.334-i*.065,f'FDR {q.fdr:.3f}',fontsize=8)
    a.set_yticks(range(3),[mod(x) for x in md.program]); a.set_ylim(2.5,-.5); a.set_xlim(float(md.ci_low.min())-.03,float(md.ci_high.max())+.03); a.set_xlabel('Native NTM-score β (95% CI)')
    q=md.iloc[0]; note(f,f'{q.dataset}: n={q.n}; MDD={q.n_target}; control={q.n_Control}; top {q.top_n}.',.17,.071,7)
    note(f,'Native ASD-derived NTM scores · secondary context · MDD limited to the available dataset.',.08,.045,7)
    note(f,'These analyses do not validate the final signed gene-PPR method.',.08,.019,7)
    export(f,name)


if __name__ == '__main__':
    (OUT / 'main').mkdir(parents=True, exist_ok=True)
    (OUT / 'supplementary').mkdir(parents=True, exist_ok=True)
    funcs = {'Figure1': figure1, 'Figure2': figure2, 'Figure3': figure3,
             'Figure4': figure4, 'Figure5': figure5,
             **{f'S{i}': fn for i, fn in enumerate([s1, s2, s3, s4, s5, s6, s7], 1)}}
    for name, fn in funcs.items():
        fn()

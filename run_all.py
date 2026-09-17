"""Generate Fig2-Fig6 and source tables with explicit reproduction modes."""
from analysis import *
import argparse,csv,json,platform
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
ROOT=Path(__file__).resolve().parent

def table(path,header,rows):
    with path.open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)

def save(fig,out,name,dpi):
    for ext in ['png','pdf','svg']:
        fig.savefig(out/f'{name}.{ext}',dpi=dpi,facecolor='white')
    plt.close(fig)

def time_axis(ax):
    ax.set(xlim=(1,541),ylim=(0,1.5),xticks=np.arange(1,542,60),xticklabels=np.arange(8,18))
    ax.grid(axis='x',color='.85',lw=.6);ax.set_axisbelow(True)

def scale_axis(ax):
    ax.set(xlim=(0,10),xticks=np.arange(1,10),xticklabels=LABELS,xlabel='Time scale (minute)')
    ax.grid(ls=':',color='.8',lw=.7);ax.set_axisbelow(True)

def right_axis(ax,label):
    twin=ax.secondary_yaxis('right',functions=(lambda x:x*300,lambda x:x/300))
    twin.set_ylabel(label)
    twin.set_yticks(ax.get_yticks()*300)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['legacy','clean'],default='clean')
    p.add_argument('--layout',choices=['single','double'],default='single',
                   help='Final print-size preset: 8.5 cm or 17.5 cm wide.')
    p.add_argument('--fig6-networks',choices=['nine','both'],default='nine',
                   help='Plot the 9-site area mean only, or retain both networks.')
    p.add_argument('--output',type=Path)
    p.add_argument('--dpi',type=int,default=600)
    p.add_argument('--width-cm',type=float,default=None,
                   help='Override the selected journal-layout width.')
    p.add_argument('--font-pt',type=float,default=7.0,
                   help='Final printed font size in points.')
    a=p.parse_args()
    preset_width={'single':8.5,'double':17.5}[a.layout]
    width_cm=a.width_cm or preset_width
    panel_height_cm={'single':6.4,'double':10.5}[a.layout]
    fig2_height_cm={'single':10.5,'double':14.0}[a.layout]
    out=a.output or ROOT/'outputs'/f'{a.mode}_{a.layout}'
    out.mkdir(parents=True,exist_ok=True)
    data=np.load(ROOT/'data'/'analysis_inputs.npz',allow_pickle=False)
    kt,site,dates=data['kt'],data['site'],data['dates'];groups=[data[k] for k in KEYS]
    assert kt.shape==(540,9,91)
    smoothed=np.stack([smooth_matlab(kt,s) for s in SPANS])
    fits=fit_correlations(smoothed,site,groups)
    stored=json.loads((ROOT/'data'/'legacy_fits.json').read_text())
    saved_coeff=np.array([[f['a'],f['b']] for f in stored])
    new_coeff=np.array([f['coefficients'] for f in fits])
    # clean changes only stale array behavior; stored model is held fixed.
    coeff=saved_coeff
    legacy=a.mode=='legacy'
    means,errors,effective,_=std_groups(smoothed,groups,legacy)
    fig6=fig6_values(smoothed,site,groups,coeff,legacy)
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
        'font.size':a.font_pt,'axes.labelsize':a.font_pt,'axes.titlesize':a.font_pt,
        'xtick.labelsize':a.font_pt,'ytick.labelsize':a.font_pt,
        'legend.fontsize':max(6.0,a.font_pt-.5),'lines.linewidth':.8,
        'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42,'ps.fonttype':42,
        'svg.fonttype':'none','axes.linewidth':.6})
    size=(width_cm/2.54,panel_height_cm/2.54)
    fig2_size=(width_cm/2.54,fig2_height_cm/2.54)
    selected=[groups[2][9],groups[1][5],groups[0][9]]
    fig,axes=plt.subplots(3,1,figsize=fig2_size,sharex=True,layout='constrained')
    for i,(ax,d) in enumerate(zip(axes,selected)):
        ax.plot(np.arange(1,541),kt[:,0,d],color='k',lw=.65)
        time_axis(ax);ax.set_ylabel('Clear-sky index')
        ax.text(.025,.86,f'({chr(97+i)}) {dates[d]}',transform=ax.transAxes)
    axes[-1].set_xlabel('Local time');save(fig,out,'Fig2',a.dpi)
    table(out/'Fig2.csv',['date','sample_index_1based','clear_sky_index'],
          ((dates[d],t+1,kt[t,0,d]) for d in selected for t in range(540)))
    fig,ax=plt.subplots(figsize=size,layout='constrained')
    for s,color,label in zip([0,2,5,8],['.8','#de6100','#0073b3','#c49900'],['1 min','10 min','1 h','4 h']):
        ax.plot(np.arange(1,541),smoothed[s,:,0,selected[0]],color=color,label=label)
    time_axis(ax);ax.set(xlabel='Local time',ylabel='Clear-sky index');ax.legend(loc='upper center',ncol=4,frameon=False)
    save(fig,out,'Fig3',a.dpi)
    table(out/'Fig3.csv',['sample_index_1based','span_1','span_11','span_61','span_241'],
          ((t+1,*smoothed[[0,2,5,8],t,0,selected[0]]) for t in range(540)))
    fig,ax=plt.subplots(figsize=size,layout='constrained')
    for c in range(3):ax.errorbar(np.arange(1,10),means[c],yerr=errors[c],color=COLORS[c],label=NAMES[c],capsize=2)
    scale_axis(ax);ax.set(ylim=(0,.3),yticks=[0,.1,.2,.3],ylabel=r'$\delta K_t$')
    right_axis(ax,r'Reference scale: $300\,\delta K_t$ (W m$^{-2}$)');ax.legend(frameon=False,loc='upper right')
    save(fig,out,'Fig4',a.dpi)
    table(out/'Fig4.csv',['sky','nominal_minutes','actual_span_samples','mean_daily_sd','sd_of_daily_sd','array_day_slots'],
          ((KEYS[c],LABELS[t],SPANS[t],means[c,t],errors[c,t],effective[c]) for c in range(3) for t in range(9)))
    fig,ax=plt.subplots(figsize=size,layout='constrained')
    for c,f in enumerate(fits):
        ax.plot(f['x'].ravel(),f['rho'].ravel(),'x',ms=2.5,mew=.6,color=COLORS[c],label=NAMES[c])
        xx=np.geomspace(np.nanmin(f['x']),np.nanmax(f['x']),400)
        ax.plot(xx,f['coefficients'][0]*xx**f['coefficients'][1],color=COLORS[c],label='Fit: '+NAMES[c].lower())
    ax.set(xscale='log',xlim=(.05,300),ylim=(0,1),xlabel='Distance/time scale (m/s)',ylabel=r'Correlation coefficient $\rho$')
    ax.grid(color='.85',lw=.6);ax.legend(frameon=False,loc='lower left');save(fig,out,'Fig5',a.dpi)
    table(out/'Fig5.csv',['sky','site_i','site_j','nominal_minutes','actual_span_samples','distance_m','distance_per_time_m_s','mean_nonnegative_rho','valid_day_count'],
          ((KEYS[c],i+1,j+1,LABELS[t],SPANS[t],f['distance'][q],f['x'][t,q],f['rho'][t,q],f['valid_days'][t,q])
           for c,f in enumerate(fits) for t in range(9) for q,(i,j) in enumerate(f['pairs'])))
    table(out/'fit_coefficients.csv',['sky','python_a','python_b','python_sse','saved_matlab_a','saved_matlab_b','saved_matlab_sse'],
          ((KEYS[c],*f['coefficients'],f['sse'],*saved_coeff[c],stored[c]['sse']) for c,f in enumerate(fits)))
    fig,ax=plt.subplots(figsize=size,layout='constrained')
    networks=[(0,9,'9-site area mean','-')]
    if a.fig6_networks=='both': networks.append((1,6,'6-site subset','--'))
    for z,n,area,linestyle in networks:
        for c in range(3):
            label=NAMES[c] if a.fig6_networks=='nine' else f'{area}: {NAMES[c].lower()}'
            ax.plot(np.arange(1,10),fig6[z,c],marker='x',ms=3,mew=.7,
                    ls=linestyle,color=COLORS[c],label=label)
    scale_axis(ax);ax.set(ylim=(0,.22),yticks=np.arange(0,.221,.02),
                          ylabel=r'SD($K_{t,\mathrm{area}}-K_{t,\mathrm{center}}$)')
    right_axis(ax,r'Reference scale: $300\times$ SD($\Delta K_t$) (W m$^{-2}$)')
    ax.legend(frameon=False,loc='upper right')
    save(fig,out,'Fig6',a.dpi)
    table(out/'Fig6.csv',['site_count','sky','nominal_minutes','actual_span_samples','delta_kt','delta_g_assuming_300_W_m2'],
          ((n,KEYS[c],LABELS[t],SPANS[t],fig6[z,c,t],300*fig6[z,c,t])
           for z,n,_,_ in networks for c in range(3) for t in range(9)))
    plotted_indices=[z for z,_,_,_ in networks]
    plotted_counts=[n for _,n,_,_ in networks]
    report=dict(mode=a.mode,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,
        kt_shape=list(kt.shape),selected_dates=dates[selected].tolist(),group_days=[len(x) for x in groups],
        fig4_array_day_slots=effective,
        fig6_first_scale={str(n):fig6[z,:,0].tolist() for z,n,_,_ in networks},
        fit_max_coefficient_difference=float(np.max(np.abs(new_coeff-saved_coeff))),
        python_fit_sse=[f['sse'] for f in fits],saved_fit_sse=[f['sse'] for f in stored],
        fig6_model='Archived MATLAB power1 fits; diagonal omitted; source distance formula retained',
        layout=a.layout,font_pt=a.font_pt,fig6_networks=a.fig6_networks,
        dimensions_cm={'Fig2':[width_cm,fig2_height_cm],
                       'Fig3-Fig6':[width_cm,panel_height_cm]},dpi=a.dpi)
    (out/'run_report.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    np.savez_compressed(out/'computed_arrays.npz',fig4_mean=means,fig4_error=errors,
                        fig6=fig6[plotted_indices],fig6_site_counts=plotted_counts,
                        coefficients=new_coeff)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()

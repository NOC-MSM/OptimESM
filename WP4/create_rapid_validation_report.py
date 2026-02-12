"""
Module containing code to generate .pdf validation reports
using RAPID-MOCHA observations.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)

"""
import os
import sys
import glob
import logging
import numpy as np
import xarray as xr
from fpdf import FPDF
from fpdf.enums import XPos, YPos
from scipy import stats
import matplotlib.pyplot as plt


# === Logging === #
def initialise_logging():
    """
    Initialise logging configuration.
    """
    logging.basicConfig(
        stream=sys.stdout,
        format="|  METRIC  | %(levelname)10s | %(asctime)s | %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# === Plotting Functions === #
# Define global figure style (Nature-style):
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.top": False,
    "ytick.right": False,
})


def _plot_amoc_streamfunctions(
    ds_mdl:xr.Dataset,
    ds_obs:xr.Dataset | None = None,
    name:str = 'Model',
    outfpath:str = 'metric_amoc_streamfunctions_26N.jpg',
    ensemble:bool = False
    ) -> None:
    """
    Plot RAPID-MOCHA 26.5°N time-mean vertical overturning stream functions.

    Parameters
    ----------
    ds_mdl: xr.Dataset
        Model dataset containing RAPID-MOCHA 26.5°N vertical overturning stream functions.

    ds_obs: xr.Dataset | None, optional
        Observational dataset containing RAPID-MOCHA 26.5°N vertical overturning stream functions.

    name: str
        Name of the model.

    outfpath: str, default='metric_amoc_streamfunctions_26N.jpg'
    Path to the plot output file.

    ensemble: bool, default=False
        Plot model ensemble average and spread.

    Returns
    -------
    None
    """
    # Define colors for plotting:
    c1='#a6cee3'
    c2='#1f78b4'

    # Extract variables from data objects
    z = ds_mdl['z']
    if ensemble:
        sf_rapid = ds_mdl['sf_rapid'].mean(dim=["time", "ens"])
        sf_model = ds_mdl['sf_model'].mean(dim=["time", "ens"])
    else:
        sf_rapid = ds_mdl['sf_rapid'].mean(dim="time")
        sf_model = ds_mdl['sf_model'].mean(dim="time")

    sfmax_rapid = sf_rapid.max()
    zmax_rapid = z[sf_rapid.argmax(dim='z')]
    sfmax_model = sf_model.max()
    zmax_model = z[sf_model.argmax(dim='z')]

    model_label = f'{name.split(" ")[0]} velocities\n(AMOC={sfmax_model:.1f} Sv, z$_{{MOC}}$={zmax_model:.0f} m)'
    rapid_label = f'{name.split(" ")[0]} RAPID approx.\n(AMOC={sfmax_rapid:.1f} Sv, z$_{{MOC}}$={zmax_rapid:.0f} m)'

    # Add data to axes:
    fig = plt.figure(figsize=(6, 8))
    plt.grid(True, ls='--', lw=1, color='0.1', alpha=0.2)

    # Model Velocities:
    if ensemble:
        plt.fill_betweenx(
            -z,
            sf_model - ds_mdl['sf_model'].mean(dim='time').std(dim='ens'),
            sf_model + ds_mdl['sf_model'].mean(dim='time').std(dim='ens'),
            color=c1, alpha=0.2)
        plt.plot(sf_model, -z, color=c1, linewidth=3, label=model_label)
    else:
        plt.plot(sf_model, -z, color=c1, linewidth=3, label=model_label)

    # RAPID Approximation
    if ensemble:
        plt.fill_betweenx(
            -z,
            sf_rapid - ds_mdl['sf_rapid'].mean(dim='time').std(dim='ens'),
            sf_rapid + ds_mdl['sf_rapid'].mean(dim='time').std(dim='ens'),
            color=c2, alpha=0.2)
        plt.plot(sf_rapid, -z, color=c2, linewidth=3, label=rapid_label)
    else:
        plt.plot(sf_rapid, -z, color=c2, linewidth=3, label=rapid_label)

    # Optional - plot observational data:
    if ds_obs is not None:
        z_obs = ds_obs.depth
        sf_obs = ds_obs.stream_function_mar
        sfmax_obs = sf_obs.mean(dim="time").max(dim="depth")
        zmax_obs = z_obs[sf_obs.mean(dim="time").argmax(dim='depth')]
        obs_label = f'RAPID Obs.\n(AMOC={sfmax_obs:.1f} Sv, z$_{{MOC}}$={zmax_obs:.0f} m)'
        plt.fill_betweenx(-z_obs, sf_obs.mean(dim="time") - sf_obs.std(dim='time'), sf_obs.mean(dim="time") + sf_obs.std(dim='time'), color='k', alpha=0.2)
        plt.plot(sf_obs.mean(dim="time"), -z_obs, 'k', linewidth=3, label=obs_label)

    # Add plot annotations:
    plt.ylim([-6200, 10])
    plt.title(f'{name}:\nVertical overturning streamfunction at 26.5°N', loc='left', fontsize=13, fontweight='bold')
    plt.xlabel('$\\psi_{_{z}}$ [Sv]', fontsize=12, fontweight='bold')
    plt.ylabel('Depth [m]', fontsize=12, fontweight='bold')
    plt.legend(loc='best', fontsize=10)

    # Save plot:
    plt.savefig(outfpath, dpi=800, bbox_inches="tight")
    plt.close(fig)


def _plot_fc_timeseries(
    ds_mdl: xr.Dataset,
    ds_obs_oht: xr.Dataset | None = None,
    ds_obs_fc: xr.Dataset | None = None,
    name: str = 'Model',
    outfpath: str = 'metric_fc_timeseries_26N.jpg',
    ensemble: bool = False
    ) -> None:
    """
    Plot Florida Current volume transport and potential temperature time series.
    
    Parameters
    ----------

    ds_mdl: xr.Dataset
    Model dataset containing RAPID 26.5°N AMOC volume transport components.

    ds_obs_oht: xr.Dataset | None, optional
    Observational dataset containing Ocean Heat Transport.

    ds_obs_fc: xr.Dataset | None, optional
    Observational dataset containing Florida Current volume transport.

    name: str, default='Model'
    Name of the model.

    outfpath: str, default='metric_fc_timeseries_26N.jpg'
    Path to the plot output file.

    ensemble: bool, default=False
        Plot model ensemble average and spread.

    Returns
    -------
    None
    """
    fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=True, constrained_layout=True)

    # Panel (a): OGCM Florida Current Volume Transport.
    axs[0].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

    fc_label = f'{name}: Florida Current = {(ds_mdl['fc'].mean()).values.item():.1f} Sv'
    if ensemble:
        axs[0].fill_between(ds_mdl['time'],
                            ds_mdl['fc'].mean(dim='ens') - ds_mdl['fc'].std(dim='ens'),
                            ds_mdl['fc'].mean(dim='ens') + ds_mdl['fc'].std(dim='ens'),
                            color="dodgerblue", alpha=0.2
                            )
        axs[0].plot(ds_mdl['time'], ds_mdl['fc'].mean(dim='ens'), linewidth=3, color="dodgerblue", label=fc_label)
    else:
        axs[0].plot(ds_mdl['time'], ds_mdl['fc'], linewidth=3, color="dodgerblue", label=fc_label)

    # Add optional observational data to sub-axis:
    if ds_obs_fc is not None:
        fc_obs_label = f'Observed: Florida Current = {ds_obs_fc.florida_current_transport.mean(skipna=True).values.item():.1f} Sv'
        axs[0].plot(ds_obs_fc['time'], ds_obs_fc.florida_current_transport, linewidth=1.5, color="0.1", alpha=0.4)
        fc_transport_monthly = ds_obs_fc.florida_current_transport.resample(time="1ME").mean()
        axs[0].plot(fc_transport_monthly['time'], fc_transport_monthly, linewidth=2.5, color="0.1", label=fc_obs_label)

    axs[0].set_ylim([np.min([ds_mdl['fc'].min(), ds_obs_fc.florida_current_transport.min()]) - 2,
                     np.max([ds_mdl['fc'].max(), ds_obs_fc.florida_current_transport.max()]) + 2]
                    )
    axs[0].set_ylabel('Volume Transport [Sv]', fontsize=12, fontweight='bold')
    axs[0].set_title(f'(a) {name}:\nFlorida Current Volume Transport at 26.5°N', loc='left', fontsize=12, fontweight='bold')
    axs[0].legend(loc=8, fontsize=10, ncol=2)

    # Panel (b): OGCM Florida Current Flow-Weighted Potential Temperature.
    axs[1].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)
    # Calculate flow-weighted temperature as: (Q_fc / (rho*Cp*T_fc))
    t_fc_fwt = ((ds_mdl['q_fc'] * 1E15) / (ds_mdl.attrs["rhocp"] * ds_mdl['fc'] * 1E6))
    fc_temp_label = f'{name}: Florida Current = {(t_fc_fwt.mean()).values.item():.1f} C'
    if ensemble:
        axs[1].fill_between(ds_mdl['time'],
                            t_fc_fwt.mean(dim='ens') - t_fc_fwt.std(dim='ens'),
                            t_fc_fwt.mean(dim='ens') + t_fc_fwt.std(dim='ens'),
                            color="dodgerblue", alpha=0.2
                            )
        axs[1].plot(ds_mdl['time'], t_fc_fwt.mean(dim='ens'), linewidth=3, color="dodgerblue", label=fc_temp_label)
    else:
        axs[1].plot(ds_mdl['time'], t_fc_fwt, linewidth=3, color="dodgerblue", label=fc_temp_label)

    # Add optional observational data to sub-axis:
    if ds_obs_oht is not None:
        fc_obs_temp_label = f'Observed: Florida Current = {ds_obs_oht.T_fc_fwt.mean(skipna=True).values.item():.1f} C'
        axs[1].plot(ds_obs_oht['time'], ds_obs_oht.T_fc_fwt, linewidth=1.5, color="0.1", alpha=0.4)
        T_fc_fwt_monthly = ds_obs_oht.T_fc_fwt.resample(time="1ME").mean()
        axs[1].plot(T_fc_fwt_monthly['time'], T_fc_fwt_monthly, linewidth=2.5, color="0.1", label=fc_obs_temp_label)

    axs[1].set_xlabel('Time', fontsize=12, fontweight='bold')
    axs[1].set_ylim([np.min([t_fc_fwt.min(), ds_obs_oht.T_fc_fwt.min()]) - 1,
                     np.max([t_fc_fwt.max(), ds_obs_oht.T_fc_fwt.max()]) + 1])
    axs[1].set_ylabel('Flow-Weighted Potential Temperature [C]', fontsize=12, fontweight='bold')
    axs[1].set_title(f'(b) {name}:\nFlorida Current Flow-Weighted Potential Temperature at 26.5°N', loc='left', fontsize=12, fontweight='bold')
    axs[1].legend(loc=8, fontsize=10, ncol=2)

    # Save plot
    plt.savefig(outfpath, dpi=800, bbox_inches="tight")
    plt.close(fig)


def _plot_amoc_transport_components(
    ds_mdl: xr.Dataset,
    ds_obs_vol: xr.Dataset | None = None,
    name: str = 'Model',
    outfpath: str = 'metric_amoc_transports_26N.jpg',
    ensemble: bool = False
    ) -> None:
    """
    Plot RAPID 26.5°N AMOC volume transport component time series.
    
    Parameters
    ----------

    ds_mdl: xr.Dataset
    Model dataset containing RAPID 26.5°N AMOC volume transport components.

    ds_obs: xr.Dataset | None, optional
    Observational dataset containing RAPID 26.5°N AMOC volume transport components.

    name: str, default='Model'
    Name of the model.

    outfpath: str, default='metric_amoc_transports_26N.jpg'
    Path to the plot output file.

    ensemble: bool, default=False
    Plot model ensemble average and spread.

    Returns
    -------
    None
    """
    fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=True, constrained_layout=True)

    # Panel (a): OGCM RAPID 26.5°N AMOC volume transport components.
    axs[0].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

    fc_label = f'Florida Current = {(ds_mdl['fc'].mean()).values.item():.1f} Sv'
    ek_label = f'Ekman = {(ds_mdl['ekman'].mean()).values.item():.1f} Sv'
    umo_label = f'Upper-Mid Ocean = {(ds_mdl['umo'].mean()).values.item():.1f} Sv'
    moc_label = f'AMOC = {(ds_mdl['moc_rapid'].mean()).values.item():.1f} Sv'

    if ensemble:
        axs[0].fill_between(ds_mdl['time'], ds_mdl['fc'].mean(dim='ens') - ds_mdl['fc'].std(dim='ens'),
                            ds_mdl['fc'].mean(dim='ens') + ds_mdl['fc'].std(dim='ens'), color="dodgerblue", alpha=0.2)
        axs[0].plot(ds_mdl['time'], ds_mdl['fc'].mean(dim='ens'), linewidth=3, color="dodgerblue", label=fc_label)

        axs[0].fill_between(ds_mdl['time'], ds_mdl['ekman'].mean(dim='ens') - ds_mdl['ekman'].std(dim='ens'),
                            ds_mdl['ekman'].mean(dim='ens') + ds_mdl['ekman'].std(dim='ens'), color="coral", alpha=0.2)
        axs[0].plot(ds_mdl['time'], ds_mdl['ekman'].mean(dim='ens'), linewidth=3, color="coral", label=ek_label)

        axs[0].fill_between(ds_mdl['time'], ds_mdl['umo'].mean(dim='ens') - ds_mdl['umo'].std(dim='ens'),
                            ds_mdl['umo'].mean(dim='ens') + ds_mdl['umo'].std(dim='ens'), color="purple", alpha=0.2)
        axs[0].plot(ds_mdl['time'], ds_mdl['umo'].mean(dim='ens'), linewidth=3, color="purple", label=umo_label)

        axs[0].fill_between(ds_mdl['time'], ds_mdl['moc_rapid'].mean(dim='ens') - ds_mdl['moc_rapid'].std(dim='ens'),
                            ds_mdl['moc_rapid'].mean(dim='ens') + ds_mdl['moc_rapid'].std(dim='ens'), color="0.1", alpha=0.2)
        axs[0].plot(ds_mdl['time'], ds_mdl['moc_rapid'].mean(dim='ens'), linewidth=3, color="0.1", label=moc_label)
    else:
        axs[0].plot(ds_mdl['time'], ds_mdl['fc'], linewidth=3, color="dodgerblue", label=fc_label)
        axs[0].plot(ds_mdl['time'], ds_mdl['ekman'], linewidth=3, color="coral", label=ek_label)
        axs[0].plot(ds_mdl['time'], ds_mdl['umo'], linewidth=3, color="purple", label=umo_label)
        axs[0].plot(ds_mdl['time'], ds_mdl['moc_rapid'], linewidth=3, color="0.1", label=moc_label)

    axs[0].set_ylim([np.min([ds_mdl['fc'].min(), ds_mdl['ekman'].min(), ds_mdl['umo'].min(), ds_mdl['moc_rapid'].min()]) - 10,
                     np.max([ds_mdl['fc'].max(), ds_mdl['ekman'].max(), ds_mdl['umo'].max(), ds_mdl['moc_rapid'].max()]) + 5]
                    )
    axs[0].set_ylabel('Volume Transport [Sv]', fontsize=12, fontweight='bold')
    axs[0].set_title(f'(a) {name}:\nVertical overturning components at 26.5°N', loc='left', fontsize=12, fontweight='bold')
    axs[0].legend(loc=8, fontsize=10, ncol=2)

    # Add optional observational data to sub-axis:
    if ds_obs_vol is not None:

        # Panel (b): RAPID-MOCHA Observations.
        axs[1].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)
        if ds_obs_vol is not None:
            fc_obs_label = f'Florida Current = {ds_obs_vol.t_gs10.mean(skipna=True).values.item():.1f} Sv'
            ek_obs_label = f'Ekman = {ds_obs_vol.t_ek10.mean(skipna=True).values.item():.1f} Sv'
            umo_obs_label = f'Upper-Mid Ocean = {ds_obs_vol.t_umo10.mean(skipna=True).values.item():.1f} Sv'
            moc_obs_label = f'AMOC = {ds_obs_vol.moc_mar_hc10.mean(skipna=True).values.item():.1f} Sv'

            axs[1].plot(ds_obs_vol['time'], ds_obs_vol.t_gs10, linewidth=1.5, color="dodgerblue", alpha=0.4)
            t_gs10_monthly = ds_obs_vol.t_gs10.resample(time="1ME").mean()
            axs[1].plot(t_gs10_monthly['time'], t_gs10_monthly, linewidth=2.5, color="dodgerblue", label=fc_obs_label)

            axs[1].plot(ds_obs_vol['time'], ds_obs_vol.t_ek10, linewidth=1.5, color="coral", alpha=0.4)
            t_ek10_monthly = ds_obs_vol.t_ek10.resample(time="1ME").mean()
            axs[1].plot(t_ek10_monthly['time'], t_ek10_monthly, linewidth=2.5, color="coral", label=ek_obs_label)
    
            axs[1].plot(ds_obs_vol['time'], ds_obs_vol.t_umo10, linewidth=1.5, color="purple", alpha=0.4)
            t_umo10_monthly = ds_obs_vol.t_umo10.resample(time="1ME").mean()
            axs[1].plot(t_umo10_monthly['time'], t_umo10_monthly, linewidth=2.5, color="purple", label=umo_obs_label)

            axs[1].plot(ds_obs_vol['time'], ds_obs_vol.moc_mar_hc10, linewidth=1.5, color="0.1", alpha=0.4)
            moc_mar_hc10_monthly = ds_obs_vol.moc_mar_hc10.resample(time="1ME").mean()
            axs[1].plot(moc_mar_hc10_monthly['time'], moc_mar_hc10_monthly, linewidth=2.5, color="0.1", label=moc_obs_label)

        axs[1].set_xlabel('Time', fontsize=12, fontweight='bold')
        axs[1].set_ylim([np.min([ds_obs_vol.t_gs10.min(), ds_obs_vol.t_ek10.min(), ds_obs_vol.t_umo10.min(), ds_obs_vol.moc_mar_hc10.min()]) - 10,
                         np.max([ds_obs_vol.t_gs10.max(), ds_obs_vol.t_ek10.max(), ds_obs_vol.t_umo10.max(), ds_obs_vol.moc_mar_hc10.max()]) + 5]
                        )
        axs[1].set_ylabel('Volume Transport [Sv]', fontsize=12, fontweight='bold')
        axs[1].set_title('(b) RAPID Observations:\nVertical overturning components at 26.5°N', loc='left', fontsize=12, fontweight='bold')
        axs[1].legend(loc=8, fontsize=10, ncol=2)

    # Save plot:
    plt.savefig(outfpath, dpi=800, bbox_inches="tight")
    plt.close(fig)


def _plot_rapid_oht_components(
    ds_mdl: xr.Dataset,
    ds_obs: xr.Dataset | None = None,
    name: str = 'Model',
    outfpath: str = 'metric_oht_components_26N.jpg',
    ensemble: bool = False
    ) -> None:
    """
    Plot RAPID ocean heat transport components time series.

    Parameters
    ----------
    ds_mdl: xr.Dataset
        Model dataset containing RAPID 26.5°N ocean heat transport components.

    ds_obs: xr.Dataset | None, optional
        Observational dataset containing RAPID 26.5°N ocean heat transport components.

    name: str, default='Model'
        Name of the model.

    outfpath: str, default='metric_oht_components_26N.jpg'
        Path to the plot output file.

    ensemble: bool, default=False
        Plot model ensemble average and spread.

    Returns
    -------
    None
    """
    # Define colors for plotting:
    c1='#a6cee3'
    c2='#1f78b4'
    c3='#b2df8a'
    c4='#33a02c'

    # Add model data to sub-axis:
    fig = plt.figure(figsize=(8,10))
    fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=True, constrained_layout=True)
    axs[0].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

    q_sum = ds_mdl['q_sum_rapid']
    q_ek = ds_mdl['q_ek']
    q_fc = ds_mdl['q_fc']
    q_geoint = ds_mdl['q_geoint']
    q_eddy = ds_mdl['q_eddy']
    q_wbw = ds_mdl['q_wbw']

    q_sum_label = f'Total = {q_sum.mean().values.item():.2f} PW'
    q_ek_label = f'Ekman = {q_ek.mean().values.item():.2f} PW'
    q_fc_label = f'Florida Current = {q_fc.mean().values.item():.2f} PW'
    q_geoint_label = f'Geostrophic Interior = {q_geoint.mean().values.item():.2f} PW'
    q_wbw_label = f'WBW = {q_wbw.mean().values.item():.2f} PW'
    q_eddy_label = f'Eddies = {q_eddy.mean().values.item():.2f} PW'

    if ensemble:
        axs[0].plot(ds_mdl['time'], q_sum.mean(dim='ens'), linewidth=3, color='k', label=q_sum_label)
        axs[0].fill_between(ds_mdl['time'], q_sum.mean(dim='ens') - q_sum.std(dim='ens'),
                            q_sum.mean(dim='ens') + q_sum.std(dim='ens'), color='k', alpha=0.2)

        axs[0].plot(ds_mdl['time'], q_ek.mean(dim='ens'), linewidth=3, color=c1, label=q_ek_label)
        axs[0].fill_between(ds_mdl['time'], q_ek.mean(dim='ens') - q_ek.std(dim='ens'),
                            q_ek.mean(dim='ens') + q_ek.std(dim='ens'), color=c1, alpha=0.2)

        axs[0].plot(ds_mdl['time'], q_fc.mean(dim='ens'), linewidth=3, color=c2, label=q_fc_label)
        axs[0].fill_between(ds_mdl['time'], q_fc.mean(dim='ens') - q_fc.std(dim='ens'),
                            q_fc.mean(dim='ens') + q_fc.std(dim='ens'), color=c2, alpha=0.2)

        axs[0].plot(ds_mdl['time'], q_wbw.mean(dim='ens'), linewidth=3, color=c3, label=q_wbw_label)
        axs[0].fill_between(ds_mdl['time'], q_wbw.mean(dim='ens') - q_wbw.std(dim='ens'),
                            q_wbw.mean(dim='ens') + q_wbw.std(dim='ens'), color=c3, alpha=0.2)

        axs[0].plot(ds_mdl['time'], q_geoint.mean(dim='ens'), linewidth=3, color=c4, label=q_geoint_label)
        axs[0].fill_between(ds_mdl['time'], q_geoint.mean(dim='ens') - q_geoint.std(dim='ens'),
                            q_geoint.mean(dim='ens') + q_geoint.std(dim='ens'), color=c4, alpha=0.2)

        axs[0].plot(ds_mdl['time'], q_eddy.mean(dim='ens'), linewidth=3, color='0.5', label=q_eddy_label)
        axs[0].fill_between(ds_mdl['time'], q_eddy.mean(dim='ens') - q_eddy.std(dim='ens'),
                            q_eddy.mean(dim='ens') + q_eddy.std(dim='ens'), color='0.5', alpha=0.2)

    else:
        axs[0].plot(ds_mdl['time'], q_sum, linewidth=3, color='k', label=q_sum_label)
        axs[0].plot(ds_mdl['time'], q_ek, linewidth=3, color=c1, label=q_ek_label)
        axs[0].plot(ds_mdl['time'], q_fc, linewidth=3, color=c2, label=q_fc_label)
        axs[0].plot(ds_mdl['time'], q_wbw, linewidth=3, color=c3, label=q_wbw_label)
        axs[0].plot(ds_mdl['time'], q_geoint, linewidth=3, color=c4, label=q_geoint_label)
        axs[0].plot(ds_mdl['time'], q_eddy, linewidth=3, color='0.5', label=q_eddy_label)

    axs[0].set_ylim([np.min([q_sum.min(), q_ek.min(), q_fc.min(), q_wbw.min(), q_geoint.min(), q_eddy.min()]) - 1,
                     np.max([q_sum.max(), q_ek.max(), q_fc.max(), q_wbw.max(), q_geoint.max(), q_eddy.max()]) + 0.5]
                    )
    axs[0].set_ylabel('Ocean Heat Transport [PW]', fontsize=12, fontweight='bold')
    axs[0].set_title(f'(a) {name}:\nOcean Heat Transports relative to 0C at 26.5°N', loc='left', fontsize=12, fontweight='bold')
    axs[0].legend(loc=8, fontsize=10, ncol=2)


    # Add observational data to sub-axis:
    if ds_obs is not None:
        axs[1].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

        q_sum_obs_label = f'Total = {(ds_obs.Q_sum.mean() / 1E15).values.item():.2f} PW'
        q_ek_obs_label = f'Ekman = {(ds_obs.Q_ek.mean() / 1E15).values.item():.2f} PW'
        q_fc_obs_label = f'Florida Current = {(ds_obs.Q_fc.mean() / 1E15).values.item():.2f} PW'
        q_geoint_obs_label = f'Geostrophic Interior = {(ds_obs.Q_int.mean() / 1E15).values.item():.2f} PW'
        q_wbw_obs_label = f'WBW = {(ds_obs.Q_wedge.mean() / 1E15).values.item():.2f} PW'
        q_eddy_obs_label = f'Eddies = {(ds_obs.Q_eddy.mean() / 1E15).values.item():.2f} PW'

        axs[1].plot(ds_obs.time, ds_obs.Q_sum / 1E15, linewidth=1.5, color='k', alpha=0.4)
        Q_sum_monthly = (ds_obs.Q_sum / 1E15).resample(time="1ME").mean()
        axs[1].plot(Q_sum_monthly.time, Q_sum_monthly, linewidth=2.5, color='k', label=q_sum_obs_label)

        axs[1].plot(ds_obs.time, ds_obs.Q_ek / 1E15, linewidth=1.5, color=c1, alpha=0.4)
        Q_ek_monthly = (ds_obs.Q_ek / 1E15).resample(time="1ME").mean()
        axs[1].plot(Q_ek_monthly.time, Q_ek_monthly, linewidth=2.5, color=c1, label=q_ek_obs_label)

        axs[1].plot(ds_obs.time, ds_obs.Q_fc / 1E15, linewidth=1.5, color=c2, alpha=0.4)
        Q_fc_monthly = (ds_obs.Q_fc / 1E15).resample(time="1ME").mean()
        axs[1].plot(Q_fc_monthly.time, Q_fc_monthly, linewidth=2.5, color=c2, label=q_fc_obs_label)

        axs[1].plot(ds_obs.time, ds_obs.Q_wedge / 1E15, linewidth=1.5, color=c3, alpha=0.4)
        Q_wedge_monthly = (ds_obs.Q_wedge / 1E15).resample(time="1ME").mean()
        axs[1].plot(Q_wedge_monthly.time, Q_wedge_monthly, linewidth=2.5, color=c3, label=q_wbw_obs_label)

        axs[1].plot(ds_obs.time, ds_obs.Q_int / 1E15, linewidth=1.5, color=c4, alpha=0.4)
        Q_int_monthly = (ds_obs.Q_int / 1E15).resample(time="1ME").mean()
        axs[1].plot(Q_int_monthly.time, Q_int_monthly, linewidth=2.5, color=c4, label=q_geoint_obs_label)

        axs[1].plot(ds_obs.time, ds_obs.Q_eddy / 1E15, linewidth=1.5, color='0.3', alpha=0.4)
        Q_eddy_monthly = (ds_obs.Q_eddy / 1E15).resample(time="1ME").mean()
        axs[1].plot(Q_eddy_monthly.time, Q_eddy_monthly, linewidth=2.5, color='0.3', label=q_eddy_obs_label)

        axs[1].set_xlabel('Time', fontsize=12, fontweight='bold')
        axs[1].set_ylim([(np.min([ds_obs.Q_sum.min(), ds_obs.Q_ek.min(), ds_obs.Q_fc.min(), ds_obs.Q_wedge.min(), ds_obs.Q_int.min(), ds_obs.Q_eddy.min()]) / 1E15) - 1,
                         (np.max([ds_obs.Q_sum.max(), ds_obs.Q_ek.max(), ds_obs.Q_fc.max(), ds_obs.Q_wedge.max(), ds_obs.Q_int.max(), ds_obs.Q_eddy.max()]) / 1E15) + 0.5]
                        )
        axs[1].set_ylabel('Ocean Heat Transport [PW]', fontsize=12, fontweight='bold')
        axs[1].set_title('(b) RAPID Observations:\nOcean Heat Transports relative to 0C at 26.5°N', loc='left', fontsize=12, fontweight='bold')
        axs[1].legend(loc=8, fontsize=10, ncol=2)

    # Save plot:
    plt.savefig(outfpath, dpi=800, bbox_inches="tight")
    plt.close(fig)


def _plot_rapid_oft_components(
    ds_mdl: xr.Dataset,
    ds_obs: xr.Dataset | None = None,
    name: str = 'Model',
    outfpath: str = 'metric_oft_components_26N.jpg',
    ensemble: bool = False
    ) -> None:
    """
    Plot RAPID ocean freshwater transport components time series.

    Parameters
    ----------
    ds_mdl: xr.Dataset
        Model dataset containing RAPID 26.5°N ocean freshwater transport components.

    ds_obs: xr.Dataset | None, optional
        Observational dataset containing RAPID 26.5°N ocean freshwater transport components.

    name: str, default='Model'
        Name of the model.

    outfpath: str, default='metric_oft_components_26N.jpg'
        Path to the plot output file.

    ensemble: bool, default=False
        Plot model ensemble average and spread.

    Returns
    -------
    None
    """
    # Define colors for plotting:
    c1='#a6cee3'
    c2='#1f78b4'
    c3='#b2df8a'
    c4='#33a02c'

    # Add model data to sub-axis:
    fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=True, constrained_layout=True)
    axs[0].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

    fw_sum = ds_mdl['fw_sum_rapid']
    fw_ek = ds_mdl['fw_ek']
    fw_fc = ds_mdl['fw_fc']
    fw_geoint = ds_mdl['fw_geoint']
    fw_eddy = ds_mdl['fw_eddy']
    fw_wbw = ds_mdl['fw_wbw']

    fw_sum_label = f'Total = {fw_sum.mean().values.item():.2f} Sv'
    fw_ek_label = f'Ekman = {fw_ek.mean().values.item():.2f} Sv'
    fw_fc_label = f'Florida Current = {fw_fc.mean().values.item():.2f} Sv'
    fw_geoint_label = f'Geostrophic Interior = {fw_geoint.mean().values.item():.2f} Sv'
    fw_wbw_label = f'WBW = {fw_wbw.mean().values.item():.2f} Sv'
    fw_eddy_label = f'Eddies = {fw_eddy.mean().values.item():.2f} Sv'

    if ensemble:
        axs[0].plot(ds_mdl['time'], fw_sum.mean(dim='ens'), linewidth=3, color='k', label=fw_sum_label)
        axs[0].fill_between(ds_mdl['time'], fw_sum.mean(dim='ens') - fw_sum.std(dim='ens'),
                            fw_sum.mean(dim='ens') + fw_sum.std(dim='ens'), color='k', alpha=0.2)

        axs[0].plot(ds_mdl['time'], fw_ek.mean(dim='ens'), linewidth=3, color=c1, label=fw_ek_label)
        axs[0].fill_between(ds_mdl['time'], fw_ek.mean(dim='ens') - fw_ek.std(dim='ens'),
                            fw_ek.mean(dim='ens') + fw_ek.std(dim='ens'), color=c1, alpha=0.2)

        axs[0].plot(ds_mdl['time'], fw_fc.mean(dim='ens'), linewidth=3, color=c2, label=fw_fc_label)
        axs[0].fill_between(ds_mdl['time'], fw_fc.mean(dim='ens') - fw_fc.std(dim='ens'),
                            fw_fc.mean(dim='ens') + fw_fc.std(dim='ens'), color=c2, alpha=0.2)

        axs[0].plot(ds_mdl['time'], fw_wbw.mean(dim='ens'), linewidth=3, color=c3, label=fw_wbw_label)
        axs[0].fill_between(ds_mdl['time'], fw_wbw.mean(dim='ens') - fw_wbw.std(dim='ens'),
                            fw_wbw.mean(dim='ens') + fw_wbw.std(dim='ens'), color=c3, alpha=0.2)

        axs[0].plot(ds_mdl['time'], fw_geoint.mean(dim='ens'), linewidth=3, color=c4, label=fw_geoint_label)
        axs[0].fill_between(ds_mdl['time'], fw_geoint.mean(dim='ens') - fw_geoint.std(dim='ens'),
                            fw_geoint.mean(dim='ens') + fw_geoint.std(dim='ens'), color=c4, alpha=0.2)

        axs[0].plot(ds_mdl['time'], fw_eddy.mean(dim='ens'), linewidth=3, color='0.5', label=fw_eddy_label)
        axs[0].fill_between(ds_mdl['time'], fw_eddy.mean(dim='ens') - fw_eddy.std(dim='ens'),
                            fw_eddy.mean(dim='ens') + fw_eddy.std(dim='ens'), color='0.5', alpha=0.2)

    else:
        axs[0].plot(ds_mdl['time'], fw_sum, linewidth=3, color='k', label=fw_sum_label)
        axs[0].plot(ds_mdl['time'], fw_ek, linewidth=3, color=c1, label=fw_ek_label)
        axs[0].plot(ds_mdl['time'], fw_fc, linewidth=3, color=c2, label=fw_fc_label)
        axs[0].plot(ds_mdl['time'], fw_wbw, linewidth=3, color=c3, label=fw_wbw_label)
        axs[0].plot(ds_mdl['time'], fw_geoint, linewidth=3, color=c4, label=fw_geoint_label)
        axs[0].plot(ds_mdl['time'], fw_eddy, linewidth=3, color='0.5', label=fw_eddy_label)

    axs[0].set_ylim([np.min([fw_sum.min(), fw_ek.min(), fw_fc.min(), fw_geoint.min(), fw_wbw.min(), fw_eddy.min()]) - 1.3,
                     np.max([fw_sum.max(), fw_ek.max(), fw_fc.max(), fw_geoint.max(), fw_wbw.max(), fw_eddy.max()]) + 0.2])
    axs[0].set_ylabel('Freshwater Transport [Sv]', fontsize=12, fontweight='bold')
    axs[0].set_title(f'(a) {name}:\nEquivalent Freshwater Transport at 26.5°N', loc='left', fontsize=12, fontweight='bold')
    axs[0].legend(loc=8, fontsize=12, ncol=2)


    # Add observational data to sub-axis:
    if ds_obs is not None:
        axs[1].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

        fw_obs_label = f'Total = {(ds_obs.frwa_trans.mean()).values.item():.2f} Sv'

        axs[1].plot(ds_obs.time, ds_obs.frwa_trans, linewidth=1.5, color='k', alpha=0.4)
        fw_obs_monthly = (ds_obs.frwa_trans).resample(time="1ME").mean()
        axs[1].plot(fw_obs_monthly.time, fw_obs_monthly, linewidth=2.5, color='k', label=fw_obs_label)

        axs[1].set_xlabel('Time', fontsize=12, fontweight='bold')
        axs[1].set_ylim([-2.8, 0])
        axs[1].set_ylabel('Freshwater Transport [Sv]', fontsize=12, fontweight='bold')
        axs[1].set_title(f'(b) {name}:\nEquivalent Freshwater Transport at 26.5°N', loc='left', fontsize=12, fontweight='bold')
        axs[1].legend(loc=8, fontsize=12, ncol=2)

    # Save plot:
    plt.savefig(outfpath, dpi=800, bbox_inches="tight")
    plt.close(fig)


def _plot_rapid_oht_geometric(
    ds_mdl: xr.Dataset,
    ds_obs: xr.Dataset,
    name: str = "Model",
    outfpath: str = 'metric_oht_geometric_26N.jpg',
    ensemble: bool = False
    ):
    """
    Plot geometric ocean heat transport components at RAPID 26.5°N.

    Parameters
    ----------
    ds_mdl : xr.Dataset
        Model dataset.
    ds_obs : xr.Dataset
        Observational dataset.
    name : str
        Name for the plot.
    outfpath : str
        Output file path for the plot.

    ensemble: bool, default=False
        Plot model ensemble average and spread.

    Returns
    -------
    None
    """
    # Add model data to sub-axis (velocities):
    fig, axs = plt.subplots(3, 1, figsize=(8, 14), sharex=True, constrained_layout=True)
    axs[0].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

    q_sum_model = ds_mdl['q_sum_model']
    q_gyre_model = ds_mdl['q_gyre_model']
    q_ot_model = ds_mdl['q_ot_model']

    q_sum_model_label = f'Total = {q_sum_model.mean().values.item():.2f} PW'
    q_ot_model_label = f'Overturning = {q_ot_model.mean().values.item():.2f} PW'
    q_gyre_model_label = f'Gyre = {q_gyre_model.mean().values.item():.2f} PW'

    if ensemble:
        axs[0].plot(ds_mdl["time"], q_sum_model.mean(dim="ens"), linewidth=3, color='k', label=q_sum_model_label)
        axs[0].fill_between(ds_mdl["time"], q_sum_model.mean(dim="ens") - q_sum_model.std(dim="ens"),
                            q_sum_model.mean(dim="ens") + q_sum_model.std(dim="ens"), color='k', alpha=0.2)

        axs[0].plot(ds_mdl["time"], q_ot_model.mean(dim="ens"), linewidth=3, color="crimson", label=q_ot_model_label)
        axs[0].fill_between(ds_mdl["time"], q_ot_model.mean(dim="ens") - q_ot_model.std(dim="ens"),
                            q_ot_model.mean(dim="ens") + q_ot_model.std(dim="ens"), color="crimson", alpha=0.2)

        axs[0].plot(ds_mdl["time"], q_gyre_model.mean(dim="ens"), linewidth=3, color="dodgerblue", label=q_gyre_model_label)
        axs[0].fill_between(ds_mdl["time"], q_gyre_model.mean(dim="ens") - q_gyre_model.std(dim="ens"),
                            q_gyre_model.mean(dim="ens") + q_gyre_model.std(dim="ens"), color="dodgerblue", alpha=0.2)

    else:
        axs[0].plot(ds_mdl["time"], q_sum_model, linewidth=3, color='k', label=q_sum_model_label)
        axs[0].plot(ds_mdl["time"], q_ot_model, linewidth=3, color="crimson", label=q_ot_model_label)
        axs[0].plot(ds_mdl["time"], q_gyre_model, linewidth=3, color="dodgerblue", label=q_gyre_model_label)

    axs[0].set_ylim([np.min([q_sum_model.min(), q_ot_model.min(), q_gyre_model.min()]) - 0.6,
                     np.max([q_sum_model.max(), q_ot_model.max(), q_gyre_model.max()]) + 0.2])
    axs[0].set_ylabel('Ocean Heat Transport [PW]', fontsize=12, fontweight='bold')
    axs[0].set_title(f'(a) {name} (velocities):\nGeometric Ocean Heat Transports at 26.5°N', loc='left', fontsize=12, fontweight='bold')
    axs[0].legend(loc=8, fontsize=12, ncol=2)

    # Add model data to sub-axis (RAPID approx):
    axs[1].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

    q_sum_rapid = ds_mdl['q_sum_rapid']
    q_gyre_rapid = ds_mdl['q_gyre_rapid']
    q_ot_rapid = ds_mdl['q_ot_rapid']

    q_sum_rapid_label = f'Total = {q_sum_rapid.mean().values.item():.2f} PW'
    q_gyre_rapid_label = f'Gyre = {q_gyre_rapid.mean().values.item():.2f} PW'
    q_ot_rapid_label = f'Overturning = {q_ot_rapid.mean().values.item():.2f} PW'

    if ensemble:
        axs[1].plot(ds_mdl["time"], q_sum_rapid.mean(dim="ens"), linewidth=3, color='k', label=q_sum_rapid_label)
        axs[1].fill_between(ds_mdl["time"], q_sum_rapid.mean(dim="ens") - q_sum_rapid.std(dim="ens"),
                            q_sum_rapid.mean(dim="ens") + q_sum_rapid.std(dim="ens"), color='k', alpha=0.2)

        axs[1].plot(ds_mdl["time"], q_ot_rapid.mean(dim="ens"), linewidth=3, color="crimson", label=q_ot_rapid_label)
        axs[1].fill_between(ds_mdl["time"], q_ot_rapid.mean(dim="ens") - q_ot_rapid.std(dim="ens"),
                            q_ot_rapid.mean(dim="ens") + q_ot_rapid.std(dim="ens"), color="crimson", alpha=0.2)

        axs[1].plot(ds_mdl["time"], q_gyre_rapid.mean(dim="ens"), linewidth=3, color="dodgerblue", label=q_gyre_rapid_label)
        axs[1].fill_between(ds_mdl["time"], q_gyre_rapid.mean(dim="ens") - q_gyre_rapid.std(dim="ens"),
                            q_gyre_rapid.mean(dim="ens") + q_gyre_rapid.std(dim="ens"), color="dodgerblue", alpha=0.2)

    else:
        axs[1].plot(ds_mdl["time"], q_sum_rapid, linewidth=3, color='k', label=q_sum_rapid_label)
        axs[1].plot(ds_mdl["time"], q_ot_rapid, linewidth=3, color="crimson", label=q_ot_rapid_label)
        axs[1].plot(ds_mdl["time"], q_gyre_rapid, linewidth=3, color="dodgerblue", label=q_gyre_rapid_label)

    axs[1].set_ylim([np.min([q_sum_rapid.min(), q_ot_rapid.min(), q_gyre_rapid.min()]) - 0.6,
                     np.max([q_sum_rapid.max(), q_ot_rapid.max(), q_gyre_rapid.max()]) + 0.2])
    axs[1].set_ylabel('Ocean Heat Transport [PW]', fontsize=12, fontweight='bold')
    axs[1].set_title(f'(b) {name} (RAPID Approx.):\nGeometric Ocean Heat Transports at 26.5°N', loc='left', fontsize=12, fontweight='bold')
    axs[1].legend(loc=8, fontsize=12, ncol=2)


    # Add optional observational data to sub-axis:
    if ds_obs is not None:
        axs[2].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)
        q_sum_obs_label = f'Total = {(ds_obs.Q_sum.mean() / 1E15).values.item():.2f} PW'
        q_gyre_obs_label = f'Gyre = {(ds_obs.Q_gyre.mean() / 1E15).values.item():.2f} PW'
        q_ot_obs_label = f'Overturning = {(ds_obs.Q_ot.mean() / 1E15).values.item():.2f} PW'

        axs[2].plot(ds_obs["time"], ds_obs.Q_sum / 1E15, linewidth=1.5, color='k', alpha=0.3)
        Q_sum_monthly = (ds_obs.Q_sum / 1E15).resample(time="1ME").mean()
        axs[2].plot(Q_sum_monthly.time, Q_sum_monthly, linewidth=2.5, color='k', label=q_sum_obs_label)

        axs[2].plot(ds_obs["time"], ds_obs.Q_ot / 1E15, linewidth=1.5, color='crimson', alpha=0.3)
        Q_ot_monthly = (ds_obs.Q_ot / 1E15).resample(time="1ME").mean()
        axs[2].plot(Q_ot_monthly.time, Q_ot_monthly, linewidth=2.5, color='crimson', label=q_ot_obs_label)

        axs[2].plot(ds_obs["time"], ds_obs.Q_gyre / 1E15, linewidth=1.5, color='dodgerblue', alpha=0.3)
        Q_gyre_monthly = (ds_obs.Q_gyre / 1E15).resample(time="1ME").mean()
        axs[2].plot(Q_gyre_monthly.time, Q_gyre_monthly, linewidth=2.5, color='dodgerblue', label=q_gyre_obs_label)

        axs[2].set_xlabel('Time', fontsize=12, fontweight='bold')
        axs[2].set_ylim([np.min([ds_obs.Q_sum.min(), ds_obs.Q_ot.min(), ds_obs.Q_gyre.min()]) / 1E15 - 0.6,
                         np.max([ds_obs.Q_sum.max(), ds_obs.Q_ot.max(), ds_obs.Q_gyre.max()]) / 1E15 + 0.2])
        axs[2].set_ylabel('Ocean Heat Transport [PW]', fontsize=12, fontweight='bold')
        axs[2].set_title('(c) RAPID Observations:\nGeometric Ocean Heat Transports at 26.5°N', loc='left', fontsize=12, fontweight='bold')
        axs[2].legend(loc=8, fontsize=12, ncol=2)


    # Save plot:
    plt.savefig(outfpath, dpi=800, bbox_inches="tight")
    plt.close(fig)


def _linreg(x, y, units='PW/Sv'):
    """ Return linear regression model and plot label """
    if len(x) > 1:
        slope, intercept, r_value, p_value, std_err =  stats.linregress(x,y)
        y_model = x * slope + intercept
        label = f'{slope:.3f} {units}'
    else:
        y_model = None
        label = ''

    return y_model, label


def _plot_amoc_vs_oht(
    ds_mdl: xr.Dataset,
    ds_obs: xr.Dataset | None = None,
    name: str = 'Model',
    outfpath: str = 'metric_amoc_vs_oht.jpg',
    ensemble: bool = False
    ) -> None:
    """
    Plot AMOC volume transports vs geometric ocean heat transport components.

    Parameters
    ----------
    ds_mdl: xr.Dataset
        Model dataset containing AMOC and OHT components.
    ds_obs: xr.Dataset | None, optional
        Observational dataset containing AMOC and OHT components.
    name: str, default='Model'
        Name of the model.
    outfpath: str, default='metric_moc_vs_oht.jpg'
        Output file path for the plot.

    ensemble: bool, default=False
        Plot model ensemble average and spread.

    Returns
    -------
    _None
    """
    # Define colors for plotting:
    c1='#a6cee3'
    c2='#1f78b4'

    fig, axs = plt.subplots(2, 1, figsize=(8, 16))

    # Add model data to axis (RAPID approximation):
    axs[0].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

    if ensemble:
        moc_rapid = ds_mdl['moc_rapid'].mean(dim='ens')
        q_sum_rapid = ds_mdl['q_sum_rapid'].mean(dim='ens')
        q_gyre_rapid = ds_mdl['q_gyre_rapid'].mean(dim='ens')
        q_ot_rapid = ds_mdl['q_ot_rapid'].mean(dim='ens')
    else:
        moc_rapid = ds_mdl['moc_rapid']
        q_sum_rapid = ds_mdl['q_sum_rapid']
        q_gyre_rapid = ds_mdl['q_gyre_rapid']
        q_ot_rapid = ds_mdl['q_ot_rapid']

    q_sum_rapid_lin, q_sum_rapid_label = _linreg(moc_rapid, q_sum_rapid)
    q_gyre_rapid_lin, q_gyre_rapid_label = _linreg(moc_rapid, q_gyre_rapid)
    q_ot_rapid_lin, q_ot_rapid_label = _linreg(moc_rapid, q_ot_rapid)

    axs[0].plot(moc_rapid, q_sum_rapid,'x', color='k', label=f'Total = {q_sum_rapid_label}')
    axs[0].plot(moc_rapid, q_ot_rapid,'x', color=c1, label=f'Overturning = {q_ot_rapid_label}')
    axs[0].plot(moc_rapid, q_gyre_rapid,'x', color=c2, label=f'Gyre = {q_gyre_rapid_label}')

    if q_sum_rapid_lin is not None:
        axs[0].plot(moc_rapid, q_sum_rapid_lin,'-', color='k')
        axs[0].plot(moc_rapid, q_ot_rapid_lin,'-', color=c1)
        axs[0].plot(moc_rapid, q_gyre_rapid_lin,'-', color=c2)

    axs[0].set_xlabel('AMOC [Sv]', fontsize=12, fontweight='bold')
    axs[0].set_ylabel('Ocean Heat transport [PW]', fontsize=12, fontweight='bold')
    axs[0].set_title(f'\n(a) {name}:\nAMOC vs OHT (RAPID approximation)', loc='left', fontsize=12, fontweight='bold')
    axs[0].legend(loc='best', fontsize=10)

    # Add optional obs data to axis:
    if ds_obs is not None:
        axs[1].grid(True, lw=0.5, ls='--', color='k', alpha=0.3)

        q_sum_obs = ds_obs.Q_sum / 1E15
        q_gyre_obs = ds_obs.Q_gyre / 1E15
        q_ot_obs = ds_obs.Q_ot / 1E15
        moc_obs = ds_obs.maxmoc

        q_sum_obs_lin, q_sum_obs_label = _linreg(moc_obs, q_sum_obs)
        q_gyre_obs_lin, q_gyre_obs_label = _linreg(moc_obs, q_gyre_obs)
        q_ot_obs_lin, q_ot_obs_label = _linreg(moc_obs, q_ot_obs)

        axs[1].plot(moc_obs, q_sum_obs,'x', color='k', label=f'Total = {q_sum_obs_label}')
        axs[1].plot(moc_obs, q_ot_obs,'x', color=c1, label=f'Overturning = {q_ot_obs_label}')
        axs[1].plot(moc_obs, q_gyre_obs,'x', color=c2, label=f'Gyre = {q_gyre_obs_label}')

        if q_sum_obs_lin is not None:
            axs[1].plot(moc_obs, q_sum_obs_lin,'-', color='k')
            axs[1].plot(moc_obs, q_ot_obs_lin,'-', color=c1)
            axs[1].plot(moc_obs, q_gyre_obs_lin,'-', color=c2)

        axs[1].set_xlabel('AMOC [Sv]', fontsize=12, fontweight='bold')
        axs[1].set_ylabel('Ocean Heat transport [PW]', fontsize=12, fontweight='bold')
        axs[1].set_title('\n(b) RAPID Observations:\nAMOC vs OHT', loc='left', fontsize=12, fontweight='bold')
        axs[1].legend(loc='best', fontsize=10)

    # Save plot:
    plt.savefig(outfpath, bbox_inches='tight', dpi=800)
    plt.close(fig)


# === Model & Observation Processing === #
def _load_rapid_obs():
    """
    Load RAPID-MOCHA observational data.

    Returns
    -------
    tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset]
        Tuple containing the in-memory RAPID-MOCHA observational datasets.
    """
    # Open Zarr stores as xr.Datasets from JASMIN Object Store:
    moc_transports = xr.open_zarr("https://noc-msm-o.s3-ext.jc.rl.ac.uk/ocean-obs/RAPID/moc_transports_v2023.1", consolidated=True)
    moc_vertical = xr.open_zarr("https://noc-msm-o.s3-ext.jc.rl.ac.uk/ocean-obs/RAPID/moc_vertical_v2023.1", consolidated=True)
    meridional_transports = xr.open_zarr("https://noc-msm-o.s3-ext.jc.rl.ac.uk/ocean-obs/RAPID/meridional_transports_v2023.1", consolidated=True)

    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    mocha_mht = xr.open_zarr("https://noc-msm-o.s3-ext.jc.rl.ac.uk/ocean-obs/RAPID/mocha_mht_data_ERA5_v2020", consolidated=True, decode_times=time_coder)
    mocha_mht['time'] = mocha_mht['time'].astype('datetime64[ns]')

    return moc_transports, moc_vertical, mocha_mht, meridional_transports


def _load_fc_obs():
    """
    Load Florida Current observational data.

    Returns
    -------
    xr.Dataset
        In-memory Florida Current observational dataset.
    """
    # Open Zarr stores as xr.Datasets from JASMIN Object Store:
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    fc_data = xr.open_zarr("https://noc-msm-o.s3-ext.jc.rl.ac.uk/ocean-obs/RAPID/fc_transport_adj_v3", consolidated=True)

    return fc_data


def _process_data(
    model_ds: xr.Dataset,
    obs_fc: xr.Dataset,
    obs_vol: xr.Dataset,
    obs_sf: xr.Dataset,
    obs_oht: xr.Dataset,
    obs_oft: xr.Dataset
    ):
    """
    Subset & load the model and observation datasets
    sharing commmon time window.

    Parameters
    ----------
    model_ds : xr.Dataset
        Model RAPID diagnostics dataset to subset.
    obs_fc : xr.Dataset
        Florida Current observation dataset.
    obs_vol : xr.Dataset
        RAPID-MOCHA volume transport dataset.
    obs_sf : xr.Dataset
        RAPID-MOCHA stream functions dataset.
    obs_oht : xr.Dataset
        RAPID-MOCHA ocean heat transport dataset.
    obs_oft : xr.Dataset
        RAPID-MOCHA ocean freshwater transport dataset.

    Returns
    -------
    tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset]
        Time subsetted model and observation datasets loaded into memory.
    """
    # Subset shared time window for model and Florida Current observations:
    fc_trange = slice(
        max([model_ds['time'].min(), obs_fc['time'].min()]),
        min([model_ds['time'].max(), obs_fc['time'].max()])
    )
    # Select only model Florida Current transport & fw-temperature variables:
    model_fc = model_ds[["fc", "q_fc"]].sel(time=fc_trange).load()
    obs_fc = obs_fc.sel(time=fc_trange).load()

    # Subset shared time window for model and RAPID-MOCHA observations:
    rapid_trange = slice(
        max([model_ds['time'].min(), obs_vol['time'].min(), obs_sf['time'].min(),
             obs_oht['time'].min(), obs_oft['time'].min()]),
        min([model_ds['time'].max(), obs_vol['time'].max(), obs_sf['time'].max(),
             obs_oht['time'].max(), obs_oft['time'].max()])
    )
    model_rapid = model_ds.sel(time=rapid_trange).load()
    obs_vol = obs_vol.sel(time=rapid_trange).load()
    obs_sf = obs_sf.sel(time=rapid_trange).load()
    obs_oht = obs_oht.sel(time=rapid_trange).load()
    obs_oft = obs_oft.sel(time=rapid_trange).load()

    return model_fc, obs_fc, model_rapid, obs_vol, obs_sf, obs_oht, obs_oft


# === Generate Validation Report === #
class ValidationReport(FPDF):
    def header(self):
        pass

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}", 0, new_x=XPos.RIGHT, new_y=YPos.TOP, align="C")


def make_validation_report(model_ds, model_name, output_pdf, ensemble=False):
    # -- Pre-Process Model & Observational Data -- #
    logging.info("In Progress: Reading RAPID 26.5N Observations")
    obs_fc = _load_fc_obs()
    obs_vol, obs_sf, obs_oht, obs_oft = _load_rapid_obs()
    logging.info("Completed: Read RAPID 26.5N Observations")

    logging.info("In Progress: Processing Model & Observation Diagnostics")
    model_fc, obs_fc, model_rapid, obs_vol, obs_sf, obs_oht, obs_oft = _process_data(model_ds, obs_fc, obs_vol, obs_sf, obs_oht, obs_oft)
    logging.info("Completed: Processed Model & Observation Diagnostics")

    if ensemble:
        dims = ["time", "ens"]
    else:
        dims = ["time"]

    pdf = ValidationReport()
    pdf.set_auto_page_break(auto=True, margin=15)

    # --- Introduction (page 0) --- #
    logging.info("In Progress: Preparing Introduction, Methodology & Summary Pages")
    pdf.add_page()
    logo_file = "./assets/METRIC_logo_v1.jpg"
    if os.path.exists(logo_file):
        pdf.image(logo_file, x=10, w=30)
    pdf.ln(1)
    pdf.set_font("Helvetica", "B", 16)
    pdf.multi_cell(0, 10,
        "Meridional Overturning Circulation Diagnostic:\n"
        "RAPID 26.5°N Validation Report",
        align="C"
    )
    pdf.ln(10)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Description", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 10)
    pdf.multi_cell(0, 6,
        "Validation report comparing Atlantic Meridional Overturning Circulation (AMOC) & "
        f"Ocean Heat and Freshwater Transport diagnostics in the {model_name} model with RAPID-MOCHA 26.5°N observations.",
        align="C",
        border="TB"
    )
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Contents", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.set_font("Helvetica", "B", 10)
    pdf.multi_cell(0, 8,
        "- Methodology\n"
        "- Summary Statistics\n"
        "- Diagnostics:",
        align="L"
    )
    diags = ["AMOC: Vertical Overtuning Streamfunctions",
             "AMOC: Volume Transport Components",
             "Florida Current: Volume Transport & Flow-Weighted Temperature",
             "OHT: Transport Components",
             "OHT: Geometric Components",
             "OFT: Transport Components",
             "AMOC-OHT: Regression",
             ]
    pdf.set_font("Helvetica", "I", 10)
    for i, diag in enumerate(diags):
        pdf.set_x(pdf.l_margin + 10)
        pdf.multi_cell(0, 8, f"{i+1}. {diag}", align="L")

    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Citation", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Please cite the associated digital object identifiers if you use METRIC for your research."
        "\n-- Code --\n"
        "Castruccio F. S., 2021: NCAR/metric: metric v0.1. doi/10.5281/zenodo.4708277"
        "\n-- Methodology--\n"
        "Danabasoglu, G., Castruccio, F. S., Small, R. J., Tomas, R., Frajka-Williams, E., and Lankhorst, M., (2021)."
        "Revisiting AMOC Transport Estimates from Observations and Models. Geophysical Research Letters.",
        align="C",
        border="TB"
    )

    # --- Introduction (page 0) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Methodology", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(5)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "METRIC enables consistent calculations of Atlantic Meridional Overturning Circulation (AMOC) estimates at the RAPID-MOCHA (26.5°N) "
        "array from observations and ocean general circulation models. To make the most appropriate comparisons, METRIC evaluates the model "
        "meridional overturning circulation using an analogous procedure to RAPID-MOCHA observations.\n"
        "\nThe RAPID array estimates the total AMOC transport across the entire Atlantic Basin at 26.5°N. The volume transport calculation "
        "involves four transport components: Florida Current (FC), western boundary wedge (WBW), mass-compensated, geostrophic mid-ocean interior (int), "
        "and the wind-driven, near-surface Ekman (EK) component (e.g., Cunningham et al., 2007; McCarthy et al., 2015).\n"
        "\nT_net = T_FC + T_WBW + T_int + T_EK\n",
        align="C",
        border="TB"
    )

    pdf.ln(4)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "User Parameters", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        f"- Florida Current (FC): Evaluated using {model_name} meridional velocities between {model_ds['v_fc'].attrs['minimum_longitude']}°E and {model_ds['v_fc'].attrs['maximum_longitude']}°E\n"
        f"- Western Boundary Wedge (WBW): Evaluated using {model_name} meridional velocities between {model_ds['wbw'].attrs['minimum_longitude']}°E and {model_ds['wbw'].attrs['maximum_longitude']}°E\n"
        f"- Geostrophic Mid-Ocean Interior (int): Evaluated using {model_name} dynamic height profiles relative to a level-of-no-motion at {model_ds.attrs['geostrophic_reference_level']} m (~{model_ds.attrs['geostrophic_reference_level']} dbar). \n"
        f"- Ekman transport: Estimated using the {model_name} zonal wind stress component. The transport is linearly distributed over the upper {model_ds.attrs['ekman_level']} m.\n"
        f"\nVertical Overturning Streamfunctions: Calculated using {model_ds.attrs['geostrophic_method']} integration following the addition of a volume (mass) balance term to ensure no net meridional flow across the RAPID 26.5°N array.",
        align="L",
        border="B"
    )

    pdf.ln(4)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Observations", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Data from the RAPID AMOC observing project is funded by the Natural Environment Research Council, U.S. National Science Foundation (NSF) with support from NOAA. "
        "They are freely available from https://rapid.ac.uk/. This validation report uses the latest RAPID observations 2023.1a (Moat et al. 2025).\n"
        "The Florida Current transport velocity times-series are funded by NOAA's Global Ocean Monitoring and Observing (GOMO) program. They are freely available from the NOAA's Western Boundary Time Series project at (https://www.aoml.noaa.gov/western-boundary-time-series/).",
        align="L",
    )

    # --- Summary (page 2) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Summary Statistics", new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.ln(5)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Atlantic Meridional Overturning Circulation at RAPID 26.5°N", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.ln(5)

    # Table: AMOC
    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(40, 5, "", 1)
    pdf.cell(45, 5, f"{model_name} (velocities)", 1)
    pdf.cell(55, 5, f"{model_name} (RAPID Approx.)", 1)
    pdf.cell(30, 5, "Observations", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("Helvetica", "", 9)
    pdf.cell(40, 5, "AMOC Max. [Sv]", 1)
    if ensemble:
        pdf.cell(45, 5, f"{model_rapid['sf_model'].max(dim='z').mean(dim=dims):.1f} ± {model_rapid['sf_model'].max(dim='z').std(dim='time').mean(dim='ens'):.1f}", 1)
        pdf.cell(55, 5, f"{model_rapid['sf_rapid'].max(dim='z').mean(dim=dims):.1f} ± {model_rapid['sf_rapid'].max(dim='z').std(dim='time').mean(dim='ens'):.1f}", 1)
    else:
        pdf.cell(45, 5, f"{model_rapid['sf_model'].max(dim='z').mean(dim=dims):.1f} ± {model_rapid['sf_model'].max(dim='z').std(dim=dims):.1f}", 1)
        pdf.cell(55, 5, f"{model_rapid['sf_rapid'].max(dim='z').mean(dim=dims):.1f} ± {model_rapid['sf_rapid'].max(dim='z').std(dim=dims):.1f}", 1)
    pdf.cell(30, 5, f"{obs_sf['stream_function_mar'].max(dim='depth').mean(dim='time'):.1f} ± {obs_sf['stream_function_mar'].max(dim='depth').std(dim='time'):.1f}", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.cell(40, 5, "Depth of AMOC Max. [m]", 1)
    pdf.cell(45, 5, f"{model_rapid['z'][model_rapid['sf_model'].mean(dim=dims).argmax(dim='z')]:.0f}", 1)
    pdf.cell(55, 5, f"{model_rapid['z'][model_rapid['sf_rapid'].mean(dim=dims).argmax(dim='z')]:.0f}", 1)
    pdf.cell(30, 5, f"{obs_sf['depth'][obs_sf['stream_function_mar'].mean(dim='time').argmax(dim='depth')]:.0f}", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.ln(5)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Atlantic Meridional Overturning Circulation: Vertical Transport Components at RAPID 26.5°N", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.ln(5)

    # Table: AMOC - Volume Transport Components [Sv]
    pdf.set_font("Helvetica", "B", 9)
    pdf.cell(40, 5, "Component", 1)
    pdf.cell(60, 5, f"{model_name}", 1)
    pdf.cell(60, 5, "Observations", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("Helvetica", "", 9)
    mdl_dict = {"moc_rapid":"Total AMOC [Sv]", "fc":"Florida Current [Sv]", "ekman":"Ekman [Sv]", "umo":"Upper-Mid Ocean [Sv]"}
    obs_dict = {"moc_rapid": "moc_mar_hc10", "fc": "t_gs10", "ekman": "t_ek10", "umo": "t_umo10"}
    for comp in mdl_dict.keys():
        mdl_mean = model_rapid[comp].mean(dim=dims).item()
        if ensemble:
            mdl_std = model_rapid[comp].std(dim='time').mean(dim='ens').item()
        else:
            mdl_std = model_rapid[comp].std(dim=dims).item()
        obs_mean = obs_vol[obs_dict[comp]].mean(dim="time").item()
        obs_std = obs_vol[obs_dict[comp]].std(dim="time").item()

        pdf.cell(40, 5, mdl_dict[comp], 1)
        pdf.cell(60, 5, f"{mdl_mean:.1f} ± {mdl_std:.1f}", 1)
        pdf.cell(60, 5, f"{obs_mean:.1f} ± {obs_std:.1f}", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.ln(5)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Ocean Heat Transport Components at RAPID 26.5°N", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.ln(5)

    # Table: Ocean Heat Transport Components [PW]
    pdf.set_font("Helvetica", "B", 9)
    pdf.cell(40, 5, "Component", 1)
    pdf.cell(65, 5, f"{model_name}", 1)
    pdf.cell(65, 5, "Observations", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("Helvetica", "", 9)
    mdl_dict = {"q_sum_rapid":"Total OHT [PW]", "q_fc":"Florida Current [PW]", "q_ek":"Ekman [PW]", "q_geoint":"Geostrophic Interior [PW]", "q_wbw":"WBW [PW]", "q_eddy":"Eddies [PW]"}
    obs_dict = {"q_sum_rapid": "Q_sum", "q_fc": "Q_fc", "q_ek": "Q_ek", "q_geoint": "Q_int", "q_wbw": "Q_wedge", "q_eddy": "Q_eddy"}
    for comp in mdl_dict.keys():
        mdl_mean = model_rapid[comp].mean(dim=dims).item()
        if ensemble:
            mdl_std = model_rapid[comp].std(dim="time").mean(dim="ens").item()
        else:
            mdl_std = model_rapid[comp].std(dim=dims).item()
        obs_mean = obs_oht[obs_dict[comp]].mean(dim="time").item() / 1E15
        obs_std = obs_oht[obs_dict[comp]].std(dim="time").item() / 1E15

        pdf.cell(40, 5, mdl_dict[comp], 1)
        pdf.cell(65, 5, f"{mdl_mean:.2f} ± {mdl_std:.2f}", 1)
        pdf.cell(65, 5, f"{obs_mean:.2f} ± {obs_std:.2f}", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.ln(5)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Ocean Freshwater Transport at RAPID 26.5°N", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
    pdf.ln(5)

    # Table: Ocean Freshwater Transport [Sv]
    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(40, 5, "", 1)
    pdf.cell(45, 5, f"{model_name} (velocities)", 1)
    pdf.cell(55, 5, f"{model_name} (RAPID Approx.)", 1)
    pdf.cell(30, 5, "Observations", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("Helvetica", "", 9)
    pdf.cell(40, 5, "Total OFT [Sv]", 1)
    if ensemble:
        pdf.cell(45, 5, f"{model_rapid['fw_sum_model'].mean(dim=dims):.1f} ± {model_rapid['fw_sum_model'].std(dim='time').mean(dim='ens'):.1f}", 1)
        pdf.cell(55, 5, f"{model_rapid['fw_sum_rapid'].mean(dim=dims):.1f} ± {model_rapid['fw_sum_rapid'].std(dim='time').mean(dim='ens'):.1f}", 1)
    else:
        pdf.cell(45, 5, f"{model_rapid['fw_sum_model'].mean(dim=dims):.1f} ± {model_rapid['fw_sum_model'].std(dim=dims):.1f}", 1)
        pdf.cell(55, 5, f"{model_rapid['fw_sum_rapid'].mean(dim=dims):.1f} ± {model_rapid['fw_sum_rapid'].std(dim=dims):.1f}", 1)
    pdf.cell(30, 5, f"{obs_oft['frwa_trans'].mean(dim='time'):.1f} ± {obs_oft['frwa_trans'].std(dim='time'):.1f}", 1, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    logging.info("Completed: Prepared Introduction, Methodology & Summary Pages")

    # --- Diagnostic Page (page 3) --- #
    logging.info("In Progress: Preparing Diagnostics Plots & Pages")
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Atlantic Meridional Overturning Circulation: Vertical Overturning Streamfunctions", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Diagnostic Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Time-mean vertical overturning streamfunction of the Atlantic Meridional "
        f"Overturning Circulation (AMOC) at 26.5°N for both {model_name} and RAPID-MOCHA observations (Moat et al., 2025)."
        f"\nTwo vertical overturning stream functions are shown for {model_name}: (1) 'velocities', which is calculated "
        "by integrating the model meridional velocity field in the longitude and depth (cumulative) dimensions, and "
        "(2) 'RAPID approximation', which is calculated using an equivalent method to the procedure used in "
        "RAPID-MOCHA observations.",
        align="L",
        border="TB"
    )
    pdf.ln(5)
    plot_file = "metric_moc_streamfunctions_26N.jpg"
    _plot_amoc_streamfunctions(ds_mdl=model_rapid, ds_obs=obs_sf, outfpath=plot_file, name=model_name, ensemble=ensemble)

    if os.path.exists(plot_file):
        pdf.image(plot_file, x=40, w=110)

    # --- Diagnostic Page (page 4) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Atlantic Meridional Overturning Circulation: Volume Transport Components", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Diagnostic Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Time-series of the Atlantic Meridional Overturning Circulation (AMOC) "
        "at 26.5°N decomposed into separate volume transport components (Florida Current, surface "
        f"Ekman layer & Upper-Mid Ocean region) for both {model_name} and RAPID-MOCHA observations (Moat et al., 2025). "
        "\nThe AMOC and its components are defined as the volume transport of the vertical overturning "
        "streamfunction at 1000 m.\n"
        "\n AMOC_Total = AMOC_fc + AMOC_ek + AMOC_umo",
        align="L",
        border="TB"
    )
    pdf.ln(5)
    plot_file = "metric_moc_components_26N.jpg"
    _plot_amoc_transport_components(ds_mdl=model_rapid, ds_obs_vol=obs_vol, outfpath=plot_file, name=model_name, ensemble=ensemble)

    if os.path.exists(plot_file):
        pdf.image(plot_file, x=35, w=130)

    # --- Diagnostic Page (page 5) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Florida Current: Volume Transport & Flow-Weighted Temperature", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Diagnostic Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Time-series of the Florida Current volume transport and flow-weighted potential temperature "
        f"for both {model_name} and RAPID-MOCHA observations (Johns et al., 2023).\n"
        "The flow-weighted potential temperature of the Florida Current is calculated by dividing the "
        "area integral of the product of velocity and potential temperature by the area integral of velocity.",
        align="L",
        border="TB"
    )
    pdf.ln(5)
    plot_file = "metric_fc_timeseries_26N.jpg"
    _plot_fc_timeseries(ds_mdl=model_fc, ds_obs_oht=obs_oht, ds_obs_fc=obs_fc, outfpath=plot_file, name=model_name, ensemble=ensemble)

    if os.path.exists(plot_file):
        pdf.image(plot_file, x=40, w=120)

    # --- Diagnostic Page (page 6) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Ocean Heat Transport: Flow Components", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Diagnostic Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Time-series of the meridional Ocean Heat Transport (OHT) "
        "at 26.5°N decomposed into the meridional temperature transports "
        "of the Florida current (Q_fc), the Ekman layer (Q_ek), the western boundary wedge (Q_wbw), "
        "the zonally averaged contribution by the mid-ocean circulation (Q_int) and the mid-ocean eddy "
        f"contribution due to spatially correlated velocity and potential temperature fluctuations (Q_eddy) for both {model_name} "
        "and RAPID-MOCHA observations (Johns et al., 2011; Johns et al., 2023).\n"
        "\n Q_Total = Q_fc + Q_ek + Q_wbw + Q_int + Q_eddy",
        align="L",
        border="TB"
    )
    pdf.ln(5)
    plot_file = "metric_oht_components_26N.jpg"
    _plot_rapid_oht_components(ds_mdl=model_rapid, ds_obs=obs_oht, name=model_name, outfpath=plot_file, ensemble=ensemble)

    if os.path.exists(plot_file):
        pdf.image(plot_file, x=35, w=130)

    # --- Diagnostic Page (page 7) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Ocean Heat Transport: Geometric Components", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Diagnostic Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Time-series of the meridional Ocean Heat Transport (OHT) at 26.5°N "
        "decomposed into geometric meridional temperature transports due to the vertical overturning "
        f"and horizontal gyre circulation for both {model_name} and RAPID-MOCHA observations (Johns et al., 2023).\n"
        "\n Q_Total = Q_ov + Q_gyre",
        align="L",
        border="TB"
    )
    pdf.ln(5)
    plot_file = "metric_oht_geometric_26N.jpg"
    _plot_rapid_oht_geometric(ds_mdl=model_rapid, ds_obs=obs_oht, name=model_name, outfpath=plot_file, ensemble=ensemble)

    if os.path.exists(plot_file):
        pdf.image(plot_file, x=45, w=110)

    # --- Diagnostic Page (page 8) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Ocean Freshwater Transport: Flow Components", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Diagnostic Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        f"Time-series of the meridional Ocean Freshwater Transport (OFT) across 26.5°N in {model_name} "
        "decomposed into the meridional freshwater transports "
        "of the Florida current (FW_fc), the Ekman layer (FW_ek), the western boundary wedge (FW_wbw), "
        "the zonally averaged contribution by the mid-ocean circulation (FW_int) and the mid-ocean eddy "
        "contribution due to spatially correlated velocity and potential temperature fluctuations (FW_eddy).\n"
        "\n FW_Total = FW_fc + FW_ek + FW_wbw + FW_int + FW_eddy\n"
        "\nThe total OFT across 26.5°N is shown for RAPID-MOCHA observations (Moat et al., 2025).",
        align="L",
        border="TB"
    )
    pdf.ln(5)
    plot_file = "metric_oft_components_26N.jpg"
    _plot_rapid_oft_components(ds_mdl=model_rapid, ds_obs=obs_oft, name=model_name, outfpath=plot_file, ensemble=ensemble)

    if os.path.exists(plot_file):
        pdf.image(plot_file, x=35, w=130)

    # --- Diagnostic Page (page 9) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "Ocean Freshwater Transport: Flow Components", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(3)
    pdf.set_font("Helvetica", "BU", 11)
    pdf.cell(0, 10, "Diagnostic Summary", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Comparison of meridional Ocean Heat Transport (OHT) and subtropical AMOC strength at 26.5°N "
        f"for both {model_name} (RAPID approximation) and RAPID-MOCHA observations (Johns et al., 2023)."
        "\n The slope (PW / Sv) of the linear model OHT = [slope]*AMOC + [intercept] is shown "
        f"for both {model_name} and RAPID-MOCHA observations.",
        align="L",
        border="TB"
    )
    pdf.ln(5)
    plot_file = "metric_amoc_vs_oht_26N.jpg"
    _plot_amoc_vs_oht(ds_mdl=model_rapid, ds_obs=obs_oht, name=model_name, outfpath=plot_file, ensemble=ensemble)

    if os.path.exists(plot_file):
        pdf.image(plot_file, x=45, w=100)

    logging.info("Completed: Prepared Diagnostics Plots & Pages")

    # --- Reference Page (page 10) --- #
    pdf.add_page()
    pdf.image(logo_file, x=10, w=20)
    pdf.ln(3)
    pdf.set_font("Helvetica", "B", 14)
    pdf.cell(0, 10, "References", new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 6,
        "Castruccio F. S., 2021: NCAR/metric: metric v0.1. doi/10.5281/zenodo.4708277\n"
        "\nDanabasoglu, G., Castruccio, F. S., Small, R. J., Tomas, R., Frajka-Williams, E., and Lankhorst, M. 2021. "
        "Revisiting AMOC Transport Estimates from Observations and Models. Geophysical Research Letters.\n"
        "\nJohns WE et al. 2011. Continuous, array-based estimates of Atlantic Ocean heat transport at 26.5°N. J. Clim. 24, 2429-2449. doi:10.1175/2010JCLI3997.1\n"
        "\nJohns William E.,  Elipot Shane, Smeed David A.,  Moat Ben, King Brian,  Volkov Denis L. and Smith Ryan H. 2023. Towards two decades of Atlantic Ocean mass and heat transports at 26.5°N Phil. Trans. R. Soc. A.38120220188. https://doi.org/10.1098/rsta.2022.0188\n"
        "\nJohns W.E., Elipot S., Smeed D.A., Moat B., King B., Volkov D.L., Smith R.H.. 2023. Atlantic Meridional Overturning Circulation (AMOC) Heat Transport Time Series between April 2004 and December 2020 at 26.5°N (v.2020) [Dataset]. University of Miami Libraries. https://doi.org/10.17604/3nfq-va20\n"
        "\nMcCarthy, G. D., Smeed, D. A., Johns, W. E., Frajka-Williams, E., Moat, B. I., Rayner, D., Baringer, M. O., Meinen, C. S., Collins, J., & Bryden, H. L. 2015. Measuring the Atlantic Meridional Overturning Circulation at 26°N. Progress in Oceanography, 130, 91-111. https://doi.org/10.1016/j.pocean.2014.10.006\n"
        "\nMoat B.I.; Smeed D.A.; Rayner D.; Johns W.E.; Smith, R.; Volkov, D.; Elipot S.; Petit T.; Kajtar J.; Baringer M. O.; and Collins, J. 2025. Atlantic meridional overturning circulation observed by the RAPID-MOCHA-WBTS (RAPID-Meridional Overturning Circulation and Heatflux Array-Western Boundary Time Series) array at 26°N from 2004 to 2023 (v2023.1a) [Dataset], British Oceanographic Data Centre - Natural Environment Research Council, UK. https://doi.org/10.5285/33826d6e-801c-b0a7-e063-7086abc0b9db\n"
        "\nVolkov, D.L., Smith, R.H., Garcia, R.F. et al. Florida Current transport observations reveal four decades of steady state. Nat Commun 15, 7780 (2024). https://doi.org/10.1038/s41467-024-51879-5",
    )

    # -- Save PDF -- #
    logging.info("In Progress: Saving validation report to .pdf")
    pdf.output(output_pdf)
    logging.info(f"Completed: Saved validation report to .pdf -> {output_pdf}")


if __name__ == "__main__":
    # ===================================== Define Input Arguments ==========================================
    ensemble = True

    model_name = "UKESM1-2-LL"
    filedir = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/WP4/data/MOHC"
    filepath = f"{filedir}/UKESM1-2-LL_esm-hist_r*_1850-01-2014-12_natl_meridional_transports_at_26N.nc"
    outfilepath = f"{filedir}/UKESM1-2-LL_esm-hist_ensemble_RAPID_26N.pdf"

    # model_name = "IPSL-CM6-ESMCO2"
    # filedir = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/WP4/data/IPSL"
    # filepath = f"{filedir}/IPSL-CM6-ESMCO2_esm-hist_r*_1850-01-2014-12_natl_meridional_transports_at_26N.nc"
    # outfilepath = f"{filedir}/IPSL-CM6-ESMCO2_esm-hist_ensemble_RAPID_26N.pdf"

    # model_name = "EC-Earth3-ESM-1"
    # filedir = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/WP4/data/"
    # filepath = f"{filedir}/*/EC-Earth3-ESM-1_esm-hist_r*_1850-01-2014-12_natl_meridional_transports_at_26N.nc"
    # outfilepath = f"{filedir}/EC-Earth3-ESM-1_esm-hist_ensemble_RAPID_26N.pdf"
    # =======================================================================================================

    # -- Import METRIC RAPID 26.5N Diagnostics -- #
    initialise_logging()
    # Define cftime coder for decoding time coordinates:
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    logging.info("In Progress: Reading model RAPID 26.5N diagnostics.")
    if ensemble:
        # Combine all ensemble members into a single Dataset:
        ds_list = []
        fpaths = glob.glob(filepath)
        for mem, file in enumerate(fpaths):
            ds_mem = xr.open_dataset(file)
            ds_mem = ds_mem.assign_coords({'ens': mem+1}).expand_dims(dim='ens', axis=0)
            ds_list.append(ds_mem)
        ds_mdl = xr.concat(ds_list, dim='ens')
    else:
        # Read single Dataset:
        ds_mdl = xr.open_dataset(filepath)

    logging.info("Completed: Read model RAPID 26.5N diagnostics.")

    # -- Run Validation Report Generation -- #
    make_validation_report(
        model_ds=ds_mdl,
        model_name=model_name,
        output_pdf=outfilepath,
        ensemble=ensemble
    )

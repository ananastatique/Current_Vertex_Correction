import numpy as np
import matplotlib.pyplot as plt

class Results():
    """
    Plot the result of the computations. The __init__ method take as the only input the
    the '.npz' file produced by the computations routines.

    All the methods take as input the ax on which draw the plot.
    """

    def __init__(self, data_file):
        rho_npz_file = np.load(data_file)
        self.B = rho_npz_file['B']
        self.T = rho_npz_file['T']

        try:
            self.rho_B = rho_npz_file['rho_B']
            cmap = plt.cm.Reds
            self.color = [cmap(i) for i in np.linspace(0, cmap.N, len(self.T)+2, dtype=np.int32)][2:]
        except:
            self.rho_T = rho_npz_file['rho_T']
            cmap = plt.cm.Blues
            self.color = [cmap(i) for i in np.linspace(0, cmap.N, len(self.B)+2, dtype=np.int32)][2:]

    def set_axis(self, B_or_T, ax, scaleB=None, **kwargs):
        """Define and set the axis of the plots."""

        if B_or_T == 'B':
            if scaleB == None:
                ax.set_xlim(0, self.B[-1])
                ax.set_xlabel(r'$B\ (T)$')
            elif scaleB == 'T':
                ax.set_xlim(1, self.B[-1]/self.T[0])
                ax.set_xlabel(r'$B/T \ (T/K)$')
            elif scaleB == 'rho0':
                ax.set_xlim(0, self.B[-1]/self.rho_B[0,0,0])
                ax.set_xlabel(r'$B/\rho_{xx}^0\ (T/\mu\Omega cm)$')
        else:
            # ax.set_xlim(0,self.T[-1])
            ax.set_xlabel(r'$T \ (K)$')
        if "ylim" in kwargs:
            ax.set_ylim(kwargs['ylim'][0], kwargs['ylim'][1])
        if 'yticks' in kwargs.keys():
            if kwargs['yticks'] == 'right':
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
        if 'xlabel' in kwargs.keys():
            if kwargs['xlabel'] == 'no':
                ax.set_xlabel('')
                ax.set_xticklabels([])
        if 'ylabel' in kwargs.keys():
            if kwargs['ylabel'] == 'no':
                ax.set_ylabel('')
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(kwargs['ylabel'])
        if 'ylim' in kwargs.keys():
            ax.set_ylim(kwargs['ylim'])

    def set_plot(self, **kwargs):
        self.plot_kwargs = {}

        self.plot_kwargs['ls'] = kwargs['ls'] if 'ls' in kwargs.keys() else None
        self.plot_kwargs['lw'] = kwargs['lw'] if 'lw' in kwargs.keys() else None
        self.plot_kwargs['marker'] = kwargs['marker'] if 'marker' in kwargs.keys() else None
        self.plot_kwargs['markerfacecolor'] = kwargs['markerfacecolor'] if 'markerfacecolor' in kwargs.keys() else None
        self.plot_kwargs['ms'] = kwargs['ms'] if 'ms' in kwargs.keys() else None

    def plot_rho_ij_vs_T(self, ij, ax, legend=True, **kwargs):
        """Plot the ij-th component of the resistivity tensor (in μΩcm) against the temperature T."""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$\rho_{}\ (\mu\Omega cm)$'.format('{'+ij+'}')
        self.set_axis('T', ax, **kwargs)
        self.set_plot(**kwargs)
        rho_ind = ['xx', 'xy', 'yx', 'yy']
        for i in range(len(self.B)):
            ln, = ax.plot(self.T, self.rho_T[rho_ind.index(ij),i,:], c=self.color[i], **self.plot_kwargs)
            if legend: ln.set_label(r'$B = {}\ T$'.format(self.B[i]))
        file_name = f'rho_{ij}_vs_T_'
        return file_name
    def plot_rho_ij_vs_B(self, ij, ax, legend=True, **kwargs):
        """Plot the ij-th component of the resistivity tensor (in μΩcm) against the magnetic field B."""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$\rho_{}\ (\mu\Omega cm)$'.format('{'+ij+'}')
        self.set_axis('B', ax, **kwargs)
        rho_ind = ['xx', 'xy', 'yx', 'yy']
        self.set_plot(**kwargs)
        for i in range(len(self.T)):
            ln, = ax.plot(self.B, self.rho_B[rho_ind.index(ij),:,i], c=self.color[i], **self.plot_kwargs)
            if legend: ln.set_label(r'$T = {}\ K$'.format(self.T[i]))
        file_name = f'rho_{ij}_vs_B_'
        return file_name
    def plot_MR(self, ax, legend=True, **kwargs):
        """Plot the magnetoresistance."""

        if 'ylabel' not in kwargs.keys():
            kwargs['ylabel'] = r'$\Delta\rho/\rho_0$'
        self.set_axis('B', ax, **kwargs)
        self.set_plot(**kwargs)
        for i in range(len(self.T)):
            ln, = ax.plot(self.B, (self.rho_B[0,:,i]-self.rho_B[0,0,i])/self.rho_B[0,0,i], c=self.color[i], **self.plot_kwargs)
            if legend :
                ln.set_label(r'$T = {}\ K$'.format(self.T[i]))
        file_name = f'MR_'
        return file_name



    def plot_RH_vs_B(self, ax, legend=True, **kwargs):
        """Plot the Hall coefficient (in mm^3/C) against the magnetic field B."""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$R_H\ (10^{-10}\Omega m/T)$' #e-10 Ohm*m/T
        self.set_axis('B', ax, **kwargs)
        self.set_plot(**kwargs)
        for i in range(len(self.T)):
            ln, = ax.plot(self.B, self.rho_B[1,:,i]*100/self.B, c=self.color[i], **self.plot_kwargs)
            if legend: ln.set_label(r'$T = {}\ K$'.format(self.T[i]))
        file_name = f'RH_vs_B_'
        return file_name

    def plot_RH_vs_T(self, ax, legend=True, **kwargs):
        """Plot the Hall coefficient (in mm^3/C) against the temperature T."""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$R_H\ (10^{-10}\Omega m/T)$' #e-10 Ohm*m/T
        self.set_axis('T', ax, **kwargs)
        self.set_plot(**kwargs)
        for i in range(len(self.B)):
            ln, = ax.plot(self.T, self.rho_T[1,i,:]*100/self.B[i], c=self.color[i], **self.plot_kwargs) #  - self.rho_T[1,i,0]*100/self.B[i]
            if legend: ln.set_label(r'$B = {}\ T$'.format(self.B[i]))
        file_name = f'RH_vs_T_'
        return file_name


    def plot_deriv_scaling(self, ax, legend=True, **kwargs):
        """Scaling plot. Plot the derivative of the resistivity as a function of B/T"""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$d\rho_{xx}/dB\ (\mu\Omega cm/T)$'
        self.set_axis('B', ax, scaleB='T', **kwargs)
        self.set_plot(**kwargs)
        drho_dB = np.zeros((len(self.B), len(self.T)))
        dB = self.B[1]-self.B[0]
        for i in range(len(self.T)):
            drho_dB[:,i] = np.gradient(self.rho_B[0,:,i], dB)
            ln, = ax.plot(self.B/self.T[i], drho_dB[:,i], c=self.color[i], **self.plot_kwargs)
            if legend :
                ln.set_label(r'$T = {}\ K$'.format(self.T[i]))
        file_name = f'deriv_scaling_'
        return file_name

    def plot_drho_db_vs_b(self, ax, legend=True, **kwargs):
        """Plot the derivative of the resistivity as a function of B"""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$d\rho_{xx}/dB\ (\mu\Omega cm/T)$'
        self.set_axis('B', ax, **kwargs)
        self.set_plot(**kwargs)
        drho_dB = np.zeros((len(self.B), len(self.T)))
        dB = self.B[1]-self.B[0]
        for i in range(len(self.T)):
            drho_dB[:,i] = np.gradient(self.rho_B[0,:,i], dB)
            ln, = ax.plot(self.B, drho_dB[:,i], c=self.color[i], **self.plot_kwargs)
            if legend :
                ln.set_label(r'$T = {}\ K$'.format(self.T[i]))
        file_name = f'drho_db_vs_b_'
        return file_name


    def plot_drho_dT_vs_b(self, ax, legend=True, **kwargs):
        """Plot the derivative /dT of the resistivity as a function of B"""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$d\rho_{xx}/dT\ (\mu\Omega cm/K)$'
        self.set_axis('T', ax, **kwargs)
        self.set_plot(**kwargs)
        drho_dT = np.zeros((len(self.B), len(self.T)))
        for i in range(len(self.B)):
            drho_dT[i,:] = np.gradient(self.rho_T[0,i,:], self.T)
            ln, = ax.plot(self.T, drho_dT[i,:], c=self.color[i], **self.plot_kwargs)
            if legend :
                ln.set_label(r'$B = {}\ T$'.format(self.B[i]))
        file_name = f'drho_dT_vs_B_'
        return file_name

    def plot_Drho_ov_T_scaling(self, ax, legend=True, **kwargs):
        """Scaling plot. Plot Δρ/T against B/T."""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$\Delta\rho/T \ (\mu\Omega cm/K)$'
        self.set_axis('B', ax, scaleB='T', **kwargs)
        self.set_plot(**kwargs)
        for i in range(len(self.T)):
            ln, = ax.plot(self.B/self.T[i], (self.rho_B[0,:,i]-self.rho_B[0,0,i])/self.T[i], c=self.color[i], **self.plot_kwargs)
            if legend :
                ln.set_label(r'$T = {}\ K$'.format(self.T[i]))
        file_name = f'Drho_ov_T_scaling_'
        return file_name
    def plot_Kohler(self, ax, legend=True, **kwargs):
        """Scaling plot. Plot Δρ/ρ0 against B/ρ0."""

        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = r'$\Delta\rho/\rho_0$'
        self.set_axis('B', ax, scaleB='rho0', **kwargs)
        self.set_plot(**kwargs)
        for i in range(len(self.T)):
            ln, = ax.plot(self.B/self.rho_B[0,0,i], (self.rho_B[0,:,i]-self.rho_B[0,0,i])/self.rho_B[0,0,i],
                          c=self.color[i], **self.plot_kwargs)
            if legend : ln.set_label(r'$T = {}\ K$'.format(self.T[i]))
        file_name = f'Kohler_'
        return file_name
    def make_legend(self, fig):
        """Draw the legend on the figure."""
        # fig.subtitle(' \n ')
        ncol = 1# len(self.T)# if len(self.T)+1<=5 else len(self.T)//2+len(self.T)%2
        fig.legend(ncol=ncol, loc='outside right center', fontsize=12,
                   fancybox=False, shadow=False)

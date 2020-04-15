import os
from tqdm.notebook import tqdm
import numpy as np
from numba import njit
import plotly
import plotly.graph_objs as go


# Plotly settings
layout = go.Layout(
    scene=dict(
        aspectmode="cube",
        camera=dict(eye=dict(x=-2, y=1.5, z=1)),
        xaxis=dict(
            title="x",
            gridcolor="rgb(255, 255, 255)",
            zerolinecolor="rgb(255, 255, 255)",
            showbackground=True,
            backgroundcolor="rgb(200, 200, 230)",
        ),
        yaxis=dict(
            title="y",
            gridcolor="rgb(255, 255, 255)",
            zerolinecolor="rgb(255, 255, 255)",
            showbackground=True,
            backgroundcolor="rgb(230, 200,230)",
        ),
        zaxis=dict(
            title="u(x, y, t)",
            gridcolor="rgb(255, 255, 255)",
            zerolinecolor="rgb(255, 255, 255)",
            showbackground=True,
            backgroundcolor="rgb(230, 230,200)",
        ),
    ),
    autosize=False,
    width=800,
    height=600,
    margin=dict(r=20, b=10, l=10, t=10),
)
camera = dict(eye=dict(x=-2, y=2, z=1))


@njit
def TDMA(coeffs, F):
    """
    Tridiagonal matrix algorithm.
    :param coeffs: matrix coefficients
    :param F: right side of the equations array
    :return: solution
    """
    N = F.size
    x = np.empty(N)
    A = np.diag(coeffs, -1)
    B = np.diag(coeffs, 1)
    C = np.diag(coeffs, 0)
    alpha = np.empty(N - 1)
    alpha[0] = -B[0] / C[0]
    beta = np.empty(N - 1)
    beta[0] = F[0] / C[0]

    # Straight run
    for i in range(1, N - 1):
        alpha[i] = -B[i] / (A[i - 1] * alpha[i - 1] + C[i])
        beta[i] = (F[i] - A[i - 1] * beta[i - 1]) / (A[i - 1] * alpha[i - 1] + C[i])
    x[-1] = (F[-1] - A[-1] * beta[-1]) / (C[-1] + A[-1] * alpha[-1])

    # Return run
    for i in range(N - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]

    return x


class HeatEquationSolver2D:
    """
    2D heat equation solver
    """

    def __init__(
        self, X_START=0, X_END=2, Y_START=0, Y_END=1, T_START=0, T_END=20, N=5, M=5, J=5
    ):

        self.X_START = X_START
        self.X_END = X_END
        self.Y_START = Y_START
        self.Y_END = Y_END
        self.T_START = T_START
        self.T_END = T_END
        self.N = N
        self.M = M
        self.J = J

        self.x = np.linspace(X_START, X_END, N)
        self.y = np.linspace(Y_START, Y_END, M)
        self.t = np.linspace(T_START, T_END, J)

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dt = self.t[1] - self.t[0]

    def initialize(
        self,
        a=1,
        f=lambda x, y, t: 0,
        fi=lambda x, y: 0,
        alpha1x=0,
        alpha2x=0,
        beta2x=0,
        beta1x=0,
        mu1x=lambda y, t: 0,
        mu2x=lambda y, t: 0,
        alpha1y=0,
        alpha2y=0,
        beta2y=0,
        beta1y=0,
        mu1y=lambda x, t: 0,
        mu2y=lambda x, t: 0,
    ):
        """
        Problem parameters initialization
        """

        self.a = a  # Laplace operator coefficient
        self.f = f  # Heterogeneity - "heat source"
        self.fi = fi  # Initial condition

        self.alpha1x = alpha1x  # Coefficient for left Neumann condition for x
        self.alpha2x = alpha2x  # Coefficient for right Neumann condition for x
        self.beta1x = beta1x  # Coefficient for left Dirichlet condition for x
        self.beta2x = beta2x  # Coefficient for right Dirichlet condition for x
        self.mu1x = mu1x  # Right side of the left boundary condition for x
        self.mu2x = mu2x  # Right side of the right boundary condition for x

        self.alpha1y = alpha1y  # Coefficient for left Neumann condition for y
        self.alpha2y = alpha2y  # Coefficient for right Neumann condition for y
        self.beta1y = beta1y  # Coefficient for left Dirichlet condition for y
        self.beta2y = beta2y  # Coefficient for right Dirichlet condition for y
        self.mu1y = mu1y  # Right side of the left boundary condition for y
        self.mu2y = mu2y  # Right side of the right boundary condition for y

        self.u = np.empty((self.J, self.N, self.M))  # u(x, y, t) initialisation

        self.u[0] = [
            [self.fi(x, y) for y in self.y] for x in self.x
        ]  # Apply initial condition

    def calculate_layer(self, j):
        """
        Step from j-1 to j layer
        """

        step_x = self.dx
        step_y = self.dt
        step_t = self.dt / 2  # Divide dt by 2 for transition to the intermediate layer

        # Transition to the intermediate layer, TDMA coefficients definition
        middle_layer = [np.array([0] * self.N)] * (self.M)
        for m in range(self.M):
            coeffs = np.empty((self.N, self.N))
            # Boundary conditions
            coeffs[0][0] = self.beta1x - self.alpha1x / step_x
            coeffs[0][1] = self.alpha1x / step_x
            coeffs[-1][-2] = -self.alpha2x / step_x
            coeffs[-1][-1] = self.alpha2x / step_x + self.beta2x

            for i in range(1, self.N - 1):
                coeffs[i][i - 1] = -self.a ** 2 * step_t / step_x ** 2
                coeffs[i][i] = 2 * self.a ** 2 * step_t / step_x ** 2 + 1
                coeffs[i][i + 1] = -self.a ** 2 * step_t / step_x ** 2

            F = np.empty(self.N)
            F[0], F[-1] = (
                self.mu1x(self.y[m], self.t[j] + step_t),
                self.mu2x(self.y[m], self.t[j] + step_t),
            )
            for k in range(1, self.N - 1):
                F[k] = (
                    self.f(self.x[k], self.y[m], self.t[j] + step_t) * step_t
                    + self.u[j - 1][k][m]
                    + self.a ** 2
                    * step_t
                    / step_y ** 2
                    * (
                        self.u[j - 1][k - 1][m]
                        - 2 * self.u[j - 1][k][m]
                        + self.u[j - 1][k + 1][m]
                    )
                )

            middle_layer[m] = TDMA(coeffs, F)  # Apply TDMA

        middle_layer = np.array(middle_layer).T

        # Transition to the next layer, TDMA coefficients definition
        new_layer = [np.array([0] * self.M)] * (self.N)

        for n in range(self.N):
            coeffs = np.empty((self.M, self.M))
            # Boundary conditions
            coeffs[0][0] = self.beta1y - self.alpha1y / step_y
            coeffs[0][1] = self.alpha1y / step_y
            coeffs[-1][-2] = self.alpha2y / step_y
            coeffs[-1][-1] = self.alpha2y / step_y + self.beta2y

            for i in range(1, self.M - 1):
                coeffs[i][i - 1] = -self.a ** 2 * step_t / step_y ** 2
                coeffs[i][i] = 2 * self.a ** 2 * step_t / step_y ** 2 + 1
                coeffs[i][i + 1] = -self.a ** 2 * step_t / step_y ** 2

            F = np.empty(self.M)
            F[0], F[-1] = (
                self.mu1y(self.x[n], self.t[j] + 2 * step_t),
                self.mu2y(self.x[n], self.t[j] + 2 * step_t),
            )
            for k in range(1, self.M - 1):
                F[k] = (
                    self.f(self.x[n], self.y[k], self.t[j] + 2 * step_t) * step_t
                    + middle_layer[n][k]
                    + self.a ** 2
                    * step_t
                    / step_x ** 2
                    * (
                        middle_layer[n][k - 1]
                        - 2 * middle_layer[n][k]
                        + middle_layer[n][k + 1]
                    )
                )

            new_layer[n] = TDMA(coeffs, F)  # Apply TDMA

        new_layer = np.array(new_layer)

        return new_layer

    def solve(self):
        print("Calculating...")
        for j in tqdm(range(1, self.J)):
            new_layer = self.calculate_layer(j)
            self.u[j] = new_layer

    def plot_state(self, n=0, save=False, filename=None):
        """
        Plot the solution state at moment n% from maximum time.
        :param n: percent of maximum time when to plot the solution
        :param save: if to save on disk
        :param filename: file name to save to
        :return: plotly.graph_objs.Figure with the state of the solution
        """
        if filename is None:
            filename = f"state-{n}%"
        if not filename.endswith(".html"):
            filename += ".html"
        num = int(round(n / 100 * self.J, 0))
        if num > self.J - 1:
            num = self.J - 1

        data = [go.Surface(x=self.x, y=self.y, z=self.u[num].T)]
        fig = go.Figure(data=data, layout=layout)
        if save:
            if not os.path.exists("results"):
                os.makedirs("results")
            plotly.offline.plot(fig, filename=f"results//{filename}")
        return fig

    def plot_initial_state(self, save=False, filename=None):
        """
        Plot the initial state.
        :param save: if to save on disk
        :param filename: file name to save to
        :return: plotly.graph_objs.Figure with the initial state of the solution
        """
        if filename is None:
            filename = "initial_state"
        return self.plot_state(n=0, save=save, filename=filename)

    def show_evolution(self, save=False, filename=None):
        """
        Animate evolution of the solution.
        :param save: if to save on disk
        :param filename: file name to save to
        :return: plotly.graph_objs.Figure with the evolution of the solution
        """
        if filename is None:
            filename = "evolution"
        if not filename.endswith(".html"):
            filename += ".html"
        fig = go.Figure(
            data=go.Surface(x=self.x, y=self.y, z=self.u[0].T), layout=layout
        )

        # Scale fix
        fig.layout.scene.zaxis.range = [np.min(self.u), np.max(self.u)]
        fig.layout.scene.yaxis.range = [self.Y_START, self.Y_END]
        fig.layout.scene.xaxis.range = [self.X_START, self.X_END]
        fig.layout.coloraxis.cmin = np.min(self.u)
        fig.layout.coloraxis.cmax = np.max(self.u)
        fig.layout.scene.xaxis.autorange = False
        fig.layout.scene.yaxis.autorange = False
        fig.layout.scene.zaxis.autorange = False

        frames = []
        for j in range(self.J):
            frames.append(
                go.Frame(data=[go.Surface(x=self.x, y=self.y, z=self.u[j].T)])
            )

        fig.frames = frames

        fig.layout.updatemenus = [
            {
                "buttons": [
                    {
                        "args": [
                            None,
                            {
                                "frame": {"duration": 100, "redraw": True},
                                "fromcurrent": True,
                            },
                        ],
                        "label": "Play",
                        "method": "animate",
                    },
                    {
                        "args": [
                            [None],
                            {
                                "frame": {"duration": 0, "redraw": True},
                                "mode": "immediate",
                            },
                        ],
                        "label": "Pause",
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top",
            }
        ]
        if save:
            if not os.path.exists("results"):
                os.makedirs("results")
            plotly.offline.plot(fig, filename=f"results//{filename}")
        return fig

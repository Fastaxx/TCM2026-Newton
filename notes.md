\documentclass[aspectratio=169]{beamer}

\usetheme{Boadilla}
\usecolortheme{default}
\usefonttheme{professionalfonts}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{tikz}
\usepackage{stmaryrd}
\usepackage{appendixnumberbeamer}
\usepackage{animate}

\title[Solving Sharp-Interface Phase-Change Problems]{Solving Sharp-Interface Phase-Change Problems:\\ Newton Is All You Need!}
\author[Louis Libat]{\underline{Louis Libat} \and Can Selçuk \and Eric Chénier \and Vincent Le Chenadec}
\institute[Gustave Eiffel University]{Multiscale Modeling and Simulation Laboratory (MSME)\\
CNRS UMR 8208, Gustave Eiffel University, France}
\date{TCM Team meeting --- \today}



\begin{document}

%------------------------------------------------------------------------------
\begin{frame}
  \titlepage
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Outline}
\tableofcontents
\end{frame}

%------------------------------------------------------------------------------
\section{Motivation}

\begin{frame}{Sharp-interface multiphysics is a coupling problem}
\begin{itemize}
  \item Phase change, heat transfer, mass transfer, free boundaries
  \item Discontinuous coefficients and jump conditions across a moving interface
  \item Geometry and PDE are strongly coupled
\end{itemize}

\vspace{0.5em}
\begin{block}{}
Sharp-interface phase-change problems are not only difficult because of physics.\\
They are also difficult because of \textbf{coupling}.
\end{block}


\begin{columns}
    \column{0.5\textwidth}
    \animategraphics[loop,controls,width=\linewidth]{10}{figures/crystal-frames/frame-}{5}{199}

    \column{0.5\textwidth}
    \animategraphics[loop,controls,width=\linewidth]{10}{figures/bubble-frames/frame-}{0}{299}
\end{columns}


\end{frame}

\begin{frame}{Two-phase sharp-interface model}

Let
\[
\Omega=\Omega^-(t)\cup\Gamma(t)\cup\Omega^+(t),
\qquad
\mathbf n \text{ points from } \Omega^- \text{ to } \Omega^+ .
\]

\vspace{0.3em}
In each phase $\Omega^\pm(t)$:

\textbf{Heat}
\[
\rho^\pm c_p^\pm
\left(
\partial_t T^\pm + \mathbf u^\pm\!\cdot\!\nabla T^\pm
\right)
=
\nabla\!\cdot\!\left(k^\pm \nabla T^\pm\right)
\]

\textbf{Species}
\[
\partial_t C^\pm + \mathbf u^\pm\!\cdot\!\nabla C^\pm
=
\nabla\!\cdot\!\left(D^\pm \nabla C^\pm\right)
\]

\textbf{Flow}
\[
\rho^\pm \partial_t \mathbf u^\pm
=
\nabla\!\cdot\!\boldsymbol{\sigma}^\pm + \rho^\pm \mathbf g,
\qquad
\nabla\!\cdot\!\mathbf u^\pm = 0
\]
with
\[
\boldsymbol{\sigma}^\pm
=
-p^\pm \mathbf I
+
\mu^\pm\left(\nabla \mathbf u^\pm + \nabla \mathbf u^{\pm T}\right).
\]

\end{frame}

\begin{frame}{Interface conditions and phase change}

On the moving interface $\Gamma(t)$, let $\mathbf w_\Gamma$ be the interface velocity
and $\dot m_\Gamma$ the interfacial mass flux.

\vspace{0.3em}
\textbf{Mass balance / kinematic relation}
\[
\dot m_\Gamma
=
\rho^- \bigl(\mathbf u^- - \mathbf w_\Gamma\bigr)\!\cdot\!\mathbf n
=
\rho^+ \bigl(\mathbf u^+ - \mathbf w_\Gamma\bigr)\!\cdot\!\mathbf n
\]

\textbf{Thermal balance (Stefan condition)}
\[
\llbracket k \nabla T \cdot \mathbf n \rrbracket
=
\dot m_\Gamma L
\]

\textbf{Species balance}
\[
\llbracket D \nabla C \cdot \mathbf n \rrbracket
=
\dot m_\Gamma \, \llbracket C \rrbracket
\]

\textbf{Partition / equilibrium law}
\[
C^+ = K\,C^-
\qquad
\text{or}
\qquad
C ^\pm = C_m 
\]

\textbf{Velocity and stress jump}
\[
\llbracket \mathbf u \rrbracket = \mathbf 0,
\qquad
\llbracket \boldsymbol{\sigma}\mathbf n \rrbracket
=
\gamma \kappa \mathbf n
\]

\vspace{0.5em}
\begin{block}{}
The PDE determines the interfacial fluxes.\\
The fluxes determine the interface motion.
\end{block}

\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{What is the usual difficulty?}
\begin{columns}[T,onlytextwidth]
\column{0.48\textwidth}
\textbf{Classical approaches}
\begin{itemize}
  \item smear the interface
  \item regularize jumps
  \item decouple PDE and interface motion
  \item use weak or explicit coupling
\end{itemize}

\column{0.48\textwidth}
\textbf{What we want instead}
\begin{itemize}
  \item keep the interface \textbf{sharp}
  \item enforce jump conditions \textbf{directly}
  \item preserve local conservation
  \item solve the coupled problem \textbf{as a nonlinear system}
\end{itemize}
\end{columns}

\vspace{1em}
\begin{alertblock}{}
We do not approximate the coupling.\\
We \textbf{solve} it.
\end{alertblock}
\end{frame}

%------------------------------------------------------------------------------
\section{Part I: Fixed interfaces (Reminder)}

\begin{frame}{Part I --- Fixed interface: sharp-interface model}
\begin{itemize}
  \item Start with a fixed embedded interface $\Gamma$
  \item Two phases $\Omega^-$ and $\Omega^+$
  \item Diffusion / heat / mass transfer with jump conditions
\end{itemize}

\vspace{0.5em}
\[
-\nabla \cdot \left( \beta^\pm \nabla u^\pm \right) = f^\pm
\qquad \text{in } \Omega^\pm
\]

\[
\llbracket \alpha u \rrbracket = g_j,
\qquad
\llbracket \beta \nabla u \cdot \bm n \rrbracket = g_f
\qquad \text{on } \Gamma
\]

\vspace{0.5em}
\begin{block}{}
Before moving interfaces, the first requirement is already clear:\\
\textbf{respect the interface geometry and the jump conditions.}
\end{block}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Numerical idea: cut cells + interface unknowns}
\begin{itemize}
  \item Cartesian finite-volume discretization
  \item Geometric moments on cut cells
  \item Bulk unknowns + interface unknowns
  \item Jump conditions enforced sharply, without smearing
\end{itemize}

\vspace{1em}
\begin{equation*}
\text{bulk PDE in cut cells}
\quad + \quad
\text{interface constraints}
\quad \Longrightarrow \quad
\text{monolithic linear system}
\end{equation*}

\vspace{1em}
\begin{block}{}
The interface is not a boundary artifact.\\
It is part of the model.
\end{block}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{What fixed-interface results already show}
\begin{itemize}
  \item Sharp jumps are captured without smearing
  \item Flux continuity is enforced at the interface
  \item Good accuracy even on difficult cut-cell configurations
  \item Clear superiority over diffuse / classical approaches
\end{itemize}

\vspace{1em}
\begin{center}
  \includegraphics[width=0.75\linewidth]{example-image}
\end{center}

\begin{center}
\small Placeholder for: field / error / comparison with classical method
\end{center}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Takeaway from Part I}
\begin{block}{Key point}
Accuracy does not come from smearing the interface less.\\
It comes from \textbf{respecting the interface exactly}.
\end{block}

\vspace{1em}
\begin{center}
\Large \textbf{A sharp interface deserves a sharp discretization.}
\end{center}
\end{frame}

%------------------------------------------------------------------------------
\section{Part II: Phase change and nonlinear coupling}

\begin{frame}{From fixed interface to phase change}
\begin{itemize}
  \item With phase change, the interface is no longer given
  \item The interface velocity depends on interfacial fluxes
  \item But those fluxes depend on the PDE solution
  \item And the PDE solution depends on the interface position
\end{itemize}

\vspace{0.8em}
\[
\text{PDE solution}
\;\Longleftrightarrow\;
\text{interfacial flux}
\;\Longleftrightarrow\;
\text{interface motion}
\]

\vspace{0.8em}
\begin{alertblock}{}
Phase change is not just a diffusion problem.\\
It is a \textbf{nonlinear coupled PDE--geometry problem}.
\end{alertblock}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{The interface is an unknown}
\begin{columns}[T,onlytextwidth]
\column{0.52\textwidth}
For a Stefan-type problem,
\[
V_n \propto \llbracket q \rrbracket
\]
with
\[
q = -\beta \nabla u \cdot \bm n
\]

So:
\begin{itemize}
  \item solve PDE $\Rightarrow$ get flux
  \item get flux $\Rightarrow$ move interface
  \item move interface $\Rightarrow$ change geometry
  \item change geometry $\Rightarrow$ new PDE
\end{itemize}

\column{0.44\textwidth}
\begin{center}
\begin{tikzpicture}[>=latex, node distance=1.6cm]
\node[draw, rounded corners, align=center] (pde) {PDE\\solution};
\node[draw, rounded corners, below of=pde, align=center] (flux) {interfacial\\flux};
\node[draw, rounded corners, below of=flux, align=center] (geo) {interface\\motion};
\draw[->] (pde) -- (flux);
\draw[->] (flux) -- (geo);
\draw[->, bend right=55] (geo.east) to (pde.east);
\end{tikzpicture}
\end{center}
\end{columns}

\begin{block}{}
The interface is not a boundary condition.\\
It is an \textbf{unknown}.
\end{block}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{What we do not want}
\begin{itemize}
  \item Explicit lagging of the interface
  \item Weak coupling between physics and geometry
  \item Iterative procedures with no clear nonlinear residual
  \item Loss of conservation / robustness when the interface moves
\end{itemize}

\vspace{1em}
\begin{alertblock}{}
If the coupling is the problem, then decoupling is not the solution.
\end{alertblock}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Our strategy: solve the coupled nonlinear problem}
\begin{itemize}
  \item Unknowns:
  \[
  \text{bulk fields} \;+\; \text{interface position}
  \]
  \item Residual contains:
  \begin{itemize}
    \item PDE equations
    \item interface jump conditions
    \item Stefan / advancement condition
  \end{itemize}
  \item Solve by nonlinear iterations
\end{itemize}

\vspace{1em}
\begin{equation*}
R(U,\Gamma)=0
\end{equation*}

\begin{center}
\Large \textbf{Do not approximate the coupling. Solve it.}
\end{center}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Nonlinear iteration: the basic idea}
At iteration $k$:
\begin{enumerate}
  \item assume an interface $\Gamma^{(k)}$
  \item build cut-cell / interfacial geometry
  \item solve the PDE on this geometry
  \item evaluate interfacial fluxes
  \item build the phase-change residual
  \item correct the interface
\end{enumerate}

\vspace{0.8em}
\[
\Gamma^{(k+1)} = \Gamma^{(k)} + \delta \Gamma^{(k)}
\]

\vspace{0.8em}
\begin{block}{}
This is not ``move first, solve later''.\\
This is a \textbf{nonlinear correction loop}.
\end{block}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Why Newton?}
\begin{itemize}
  \item The difficulty is the nonlinear coupling
  \item Newton provides a systematic correction mechanism
  \item Natural framework for monolithic PDE + interface updates
  \item Same philosophy extends to richer physics
\end{itemize}

\vspace{1em}
\begin{alertblock}{}
Sharp-interface phase change is a nonlinear problem.\\
So we use a nonlinear solver.
\end{alertblock}

\vspace{0.8em}
\begin{center}
\Large \textbf{Newton is all you need.}
\end{center}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Extension: beyond pure diffusion}
The same framework extends to more coupled settings:
\begin{itemize}
  \item pure diffusion Stefan problems
  \item heat and mass transfer
  \item coupling with Stokes flow
  \item future: more general two-fluid multiphysics
\end{itemize}

\vspace{0.8em}
\begin{block}{}
Once the interface is treated as an unknown,\\
the framework becomes naturally \textbf{multiphysics-ready}.
\end{block}
\end{frame}

%------------------------------------------------------------------------------
\section{Results and code}

\begin{frame}{Illustrative results}
\begin{itemize}
  \item fixed-interface validation
  \item phase-change nonlinear convergence
  \item 1D / 2D / 3D examples
  \item diffusion and diffusion--Stokes coupling
\end{itemize}

\vspace{1em}
\begin{center}
  \includegraphics[width=0.75\linewidth]{example-image}
\end{center}

\begin{center}
\small Placeholder for: convergence plot / interface evolution / residual decay
\end{center}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}{Code perspective}
\begin{itemize}
  \item Cartesian cut-cell framework
  \item sharp-interface operators
  \item interface-aware PDE assembly
  \item nonlinear coupling between resolution and advancement
\end{itemize}

\vspace{1em}
\begin{block}{}
The code reflects the mathematical idea:\\
geometry and PDE are solved together, not separately.
\end{block}
\end{frame}

%------------------------------------------------------------------------------
\section{Conclusion}

\begin{frame}{Take-home messages}
\begin{enumerate}
  \item Sharp-interface methods already outperform classical smeared approaches for fixed interfaces
  \item For phase change, the real difficulty is the PDE--geometry coupling
  \item The interface must be treated as an unknown
  \item Nonlinear iterations provide a natural framework to solve the full problem
\end{enumerate}

\vspace{1em}
\begin{alertblock}{}
The interface is not a boundary condition.\\
It is an unknown.
\end{alertblock}

\vspace{0.6em}
\begin{center}
\Large \textbf{In Newton we trust.}
\end{center}
\end{frame}

%------------------------------------------------------------------------------
\begin{frame}
\centering
\Huge Thank you
\vspace{1em}

\Large Questions?
\end{frame}

\end{document}
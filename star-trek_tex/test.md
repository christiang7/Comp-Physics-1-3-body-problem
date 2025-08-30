===== test =====
# test
Created 2025-08-28


## Description

## Journal
 - [x] Backlog
    - [ ] 
 - [x] Doing

## Latex File


\documentclass[10pt,a4paper,listof=totoc]{scrreprt}
\usepackage[left=2.5cm,right=2.5cm,top=1cm,bottom=2cm,includeheadfoot]{geometry}
\usepackage{pslatex,palatino,avant,graphicx,color}
\usepackage[german]{babel}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}
\usepackage{float}
%\usepackage{scrhack}
\usepackage{amsfonts}
\usepackage{graphicx}
\usepackage{epstopdf}
\renewcommand{\chapterheadstartvskip}{\vspace*{.1\baselineskip}}
\usepackage[font=scriptsize]{caption}
\usepackage{subcaption}
\usepackage[format=hang]{caption}
\usepackage{setspace}
\usepackage[urlcolor=blue]{hyperref}
\hypersetup{hidelinks,colorlinks=false}
\usepackage{etoolbox}
%\usepackage{natbib}
\usepackage{multibib}
\usepackage{listings}
\usepackage{tikz}
\usepackage{comment}
\usepackage{xcolor}
\usepackage{minted}
\usemintedstyle{friendly}
\definecolor{lightgraycolor}{rgb}{0.95,0.95,0.95}

\onehalfspacing
\makeatletter
\renewcommand*{\size@chapter}{\large}
\makeatother
\setkomafont{section}{\normalfont\bfseries}                % Titel mit Normalschrift
\setkomafont{subsection}{\normalfont\bfseries}             % Titel mit Normalschrift
\setkomafont{subsubsection}{\normalfont\bfseries}          % Titel mit Normalschrift
%\setkomafont{caption}{\normalfont\large\bfseries}         % Fette Beschriftungen
%\setkomafont{captionlabel}{\normalfont\bfseries}          % Fette Beschriftungen
\captionsetup{font=normalsize}
%\captionsetup[figure]{labelfont=normalsize}

*test.tex*
```tex
\documentclass[DIN,
	% Schriftgröße
	fontsize=12pt,
a4paper]{scrartcl}
\input{general-preamble.tex}
\input{color-style.tex}
\usepackage{longtable}
\usepackage{physics}
\begin{document}
\input{./Deck.tex}
\tableofcontents
\newpage
\part{Einleitung}
Unendliche Weiten, wir schreiben das Jahr 2017 und unser Raumschiff fliegt seine Bahnen um unser Erde-Mond System. Wir werden in unserer Arbeit das Dreikörperproblem genauer untersuchen. Wir betrachten das Gravitationspotential zweier Himmelskörper wie die Erde und den Mond und lassen einen Probekörper - Raumschiff durch dieses System hindurch fliegen. Wie bereits bekannt ist das Dreikörperproblem nicht allgemein exakt analytisch lösbar. Es gibt viele Möglichkeiten das Problem näherungsweise zu lösen. Hier betrachten wir die zwei Körper Erde-Mond als massereiche Körper und ein Raumschiff, was nur eine verschwindend geringe Masse besitzt, was die Dynamik des Erde-Mond Systems kaum stört. Demnach können wir den Hamiltonoperator und daraus die Bewegungsgleichungen des Raumschiffs angeben. Es ist ein gekoppeltes Differentialgleichzugsystem, wodurch die Trajektorie des Raumschiffs stark von der Wahl der Anfangsbedingungen abhängt. Die Erde und der Mond befinden sich in einer Ebene und rotieren gemeinsam als ein starrer Körper entgegengesetzt dem Uhrzeigersinn. Wir gehen über ins mitrotierende System, wo sich die Erde bei $ (-\mu,0) $ und der Mond bei $ (1-\mu, 0) $ befindet. Demnach befindet sich der Schwerpunkt im Ursprung des Koordinatensystems. In unserer Arbeit werden die folgenden Themen näher untersuchen: Hills Region, Lagrange-Punkte, Symplektisches Intergrationsverfahren und einige interessante Flugbahnen des Raumschiff finden und darstellen.

\begin{figure}[htp!]
	\begin{center}
		\begin{tikzpicture}
		\draw[fill=blue] (0,0) node[above=1.05cm] {Erde} circle [radius=1cm];
		\draw[fill=gray] (5,0) node[above=0.25cm] {Mond} circle [radius=0.2cm];
		\draw[fill=red] (2,2) node[above] {Raumschiff} circle [radius=0.05cm];
		\end{tikzpicture}
		\caption{Das Erde-Mond System mit Raumschiff im mitrotierenden System.}
	\end{center}
	\label{pic:EMS}
\end{figure}
\newpage
\part{Grundlagen}\label{ch:grd}
Die Dynamik dieses Dreikörperproblems wird von folgenden Hamiltonian beschrieben.
\begin{align}\label{ham}
H&=\frac{(p_x + y)^2}{2}+\frac{(p_y-x)^2}{2}- \Omega(x,y)\\
\text{mit }\Omega(x,y)&=\frac{x^2+y^2}{2} +\frac{1-\mu}{r_1}+\frac{\mu}{r_2}+\frac{\mu(1-\mu)}{2}\notag\\
\text{wobei }r_1&=\sqrt{(x+\mu)^2+ y^2}, \quad r_2=\sqrt{(x-1+\mu)^2+ y^2}\notag
\end{align}
In Gleichung \eqref{ham} ist $\mu=m_M$ und $m_E=1-\mu$. Die Gesamtmasse des Erd-Mond-Systems ist somit normiert. Die Hamiltonschen Bewegungsgleichungen lauten:
\begin{align}
&\left\{
\begin{array}{l l}
\dot x=\frac{dx}{dt}=\frac{\partial H}{\partial p_x}=p_x+y\\
\dot y=\frac{dy}{dt}=\frac{\partial H}{\partial p_y}=p_y-x
\end{array}
\right.\ \label{Bew_glg_x}\\
&\left\{
\begin{array}{l l}
\dot p_x=\frac{dp_x}{dt}=-\frac{\partial H}{\partial x}=p_y-x+\Omega_x\\
\dot p_y=\frac{dp_y}{dt}=-\frac{\partial H}{\partial y}=-p_x-y+\Omega_y
\end{array}
\right.\ \label{Bew_glg_p}\\
\text{mit }\Omega_x&=\frac{\partial\Omega}{\partial x}, \quad \Omega_y=\frac{\partial\Omega}{\partial y}\notag
\end{align}
Leitet man  die Glg. \eqref{Bew_glg_x} nach $t$ ab und setzt Glg. \eqref{Bew_glg_p} ein so erhält man:
\begin{align}\label{Bew_glg}
\left\{
\begin{array}{l l}
\ddot x=\dot p_x+\dot y=p_y-x+\Omega_x+\dot y\\
\ddot y=\dot p_y-\dot x=-p_x-y+\Omega_y-\dot x
\end{array}
\right.\
\end{align}
Stellt man Glg. \eqref{Bew_glg_x} nach $p=(p_x, p_y)^T$ um und setzt diese Ausdrücke in Glg. \eqref{Bew_glg} ein, erhält man schließlich die grundlegenden Bewegungsgleichungen, die wir für unserer weiteren Untersuchungen numerisch analysieren werden.
\begin{align}\label{Bew}
\left\{
\begin{array}{l l}
\ddot x=\dot y+x-x+\Omega_x+\dot y=2\dot y+\Omega_x\\
\ddot y=-\dot x+y-y+\Omega_y-\dot x=-2\dot x+\Omega_y
\end{array}
\right.\
\end{align}
Des Weiteren ist  die Energie des Systems $E$ wie folgt definiert.
\begin{align}\label{Energie}
E=\Omega(x,y)-\frac{\dot x^2+\dot y^2}{2}
\end{align}
Bei näheren Untersuchungen stellt sich heraus, dass die Energie eine Erhaltungsgröße ist.
\begin{align}\label{Energieerhalt}
\dot E=\dot \Omega(x,y)-\ddot x \dot x-\ddot y \dot y=\Omega_x \dot x+\Omega_y \dot y-\ddot x \dot x-\ddot y \dot y
\end{align}
Aus Gleichung \eqref{Bew} erhält man $\Omega_y=\ddot y+2\dot x$ und $\Omega_x=\ddot x-2\dot y$. Setzt man dies nun in Glg. \eqref{Energieerhalt} ein erkennt man, dass die Energie erhalten ist.
\begin{align*}
\dot E&=(\ddot x-2\dot y) \dot x+(\ddot y+2\dot x) \dot y-\ddot x \dot x-\ddot y \dot y\notag\\
 &=\ddot x\dot x-2\dot y\dot x+\ddot y \dot y+2\dot x \dot y -\ddot x \dot x-\ddot y \dot y=0
\end{align*}
\newpage
Für die numerische Analyse des Differentialgleichungssystem \eqref{Bew} transformieren wir dieses zu einem Differentialgleichungssystem 1. Ordnung mit 4 Variablen:
\begin{align}\label{dfglg}
&\left\{
\begin{array}{l l l l}
\dot z_1=z_3\\
\dot z_2=z_4\\
\dot z_3=2z_4+z_1-\frac{(1-\mu)\cdot(z_1+\mu)}{r_1^3}-\frac{\mu\cdot(z_1-1+\mu)}{r_2^3}=2z_4+\Omega_{z_1}(z_1,z_2)\\
\dot z_4=-2z_3+z_2-\frac{(1-\mu)\cdot z_2}{r_1^3}-\frac{\mu \cdot z_2}{r_2^3}=-2z_3+\Omega_{z_2}(z_1,z_2)
\end{array}
\right.\ \\
\text{mit } z_1&=x \text{ und } z_2=y\notag\\
\text{wobei }r_1&=\sqrt{(z_1+\mu)^2+ z_2^2}, \quad r_2=\sqrt{(z_1-1+\mu)^2+ z_2^2}\notag
\end{align}
Um herauszufinden, wo die Lagrange-Punkte $L_1-L_5$ (Gleichgewichtspunkte, siehe Abb. \ref{pic:Lag}) liegen, erhalten wir mithilfe des Euler-Verfahren folgende Iterationsvorschrift:
\begin{align}\label{glgbed}
&\left\{
\begin{array}{l l l l}
z_1^{n+1}=z_1^n+z_3^n h\\
z_2^{n+1}=z_2^n+z_4^n h\\
z_3^{n+1}=z_3^n+\Omega_{z_1}(z_1^n,z_2^n)h\\
z_4^{n+1}=z_4^n+\Omega_{z_2}(z_1^n,z_2^n)h
\end{array}
\right.\ \notag\\
\Leftrightarrow \vec{z}^{\;n+1}&=\vec{z}^{\;n}+\vec {f} \left(z_1^n,z_2^n,z_3^n,z_4^n\right)\cdot h\\
\text{mit } h&=t^{n+1}-t^{n} \text{ und }\vec{z}^{\;n}=\left(z_1^n,z_2^n,z_3^n,z_4^n\right)^T\notag
\end{align}
Die Gleichgewichtsbedingung für Gleichung \eqref{glgbed} lautet $\vec{z}^{\;n+1}=\vec{z}^{\;n}$ und ist erfüllt, wenn $\vec {f} (z_1^n,z_2^n,z_3^n,z_4^n)=\vec 0$ gilt. Somit sind dort $z_3^n=0$ und $z_4^n=0 \; \forall \; n \in \mathbb{Z}^+$. Für $L_1, L_2$ und $L_3$ ist $z_2=0$ und somit $\Omega_{z_2}(z_1,z_2=0)=0$. Damit vereinfacht sich die Findung dieser Lagrangepukte zur einfachen Bestimmung der Nullstelle von $\Omega_{z_1}(z_1,z_2=0)=0$. Diese bestimmen wir mithilfe des Bisektionsverfahren (Glg. \eqref{gleichung bisec}, S. \pageref{gleichung bisec}).
\newpage
\part{Auswertung}
\section{Lagrangepunkte und Hill's Region} \label{kap_2}
Die Werte für die Koordinaten und die Potentiale an den Lagrangepunken sind in der Tabelle \ref{tab:Lag} aufgeführt.\\
\begin{figure}[htp!]
	\begin{center}
		\begin{tikzpicture}[scale=0.9]
		\draw[fill=gray] (5,0) circle [radius=0.1cm];
		\draw (0,0) circle [radius=1cm];
		%\filldraw[fill=gray!20,draw=gray!50!black] (0,0) -- (2cm,0cm) arc [start angle=0, end angle=60, radius=2cm] -- cycle;
		%\filldraw[fill=gray!20,draw=gray!50!black] (0,0) -- (2cm,0cm) arc [start angle=0, end angle=-60, radius=2cm] -- cycle;
		\draw[fill=red] (2.5,4.33) node[above] {$ L_4 $} circle [radius=0.05cm];
		\draw[fill=red] (2.5,-4.33) node[above] {$ L_5 $} circle [radius=0.05cm];
		\draw[fill=red] (-5,0) node[left=0.25cm,above=0.05cm] {$ L_3 $} circle  [radius=0.05cm];
		\draw[fill=red] (4.7,0) node[left=0.1cm,above] {$ L_1 $} circle [radius=0.05cm];
		\draw[fill=red] (5.3,0) node[right=0.1cm,above] {$ L_2 $} circle [radius=0.05cm];
		\draw[fill=blue] (0,0) circle [radius=1cm];
		\draw (0,0) circle [radius=5cm];
		\draw[dotted] (0,0) -- (2.5,4.33);
		\draw[dotted] (0,0) -- (2.5,-4.33);
		\draw[dotted] (0,0) -- (-5.5,0);
		\draw[dotted] (0,0) -- (5.5,0);
		\end{tikzpicture}
	\end{center}
	\caption{Die Lagrange Punkte im Gravitationsfeld des Erde-Mond Systems.}
	\label{pic:Lag}
\end{figure}\\

\begin{comment}
\begin{figure}[h]
	\centering
		\includegraphics[width=0.65\textwidth]{Konstellation.png}
		\caption{Lage der Lagrangepunkte}
	\label{fig:konst}
\end{figure}\\
\end{comment}
\begin{longtable}{|l|l|l|l|c|c|c|}
  \hline
	$i$	& $L_i=(x_i,y_i)$ 	& $\Omega_i=\Omega(x_i,y_i)$	 \\ \hline
	1		& (0.716,0) 	& 1.734\\ \hline
  2   &  (1.228,0)  & 1.701 \\ \hline
	3   &  (-1.021,0) & 1.549  \\ \hline
	 \caption{Tabelle mit Koordinaten der Lagrangepunkte $L_i$ (Werte mit einer Genauigkeit von 0.001) und der dort existenten Potentiale $\Omega(L_i)$.}
	\label{tab:Lag}
\end{longtable}
\newpage
Die Menge aller Punkte $(x,y)$ für die $\Omega(x,y)\geq E$ definiert ein Gebiet namens Hill's Region, welches für ein Raumfahrzeug mit der Energie $E$ erreichbar ist. In der Abbildung \ref{fig:Hill} sind die Begrenzungen von 6 solcher Regionen dargestellt. \\
\begin{figure}[h]
\input{Hills.tex}
	\centering
		\caption{Potential $\Omega$ um $L_2$ mit 6 Begrenzungen der jeweiligen Hill's Region.}
	\label{fig:Hill}
\end{figure}\\
Anhand der Abb. \ref{fig:Hill} ist ersichtlich, dass sich für Energien $E\leq\Omega_2$ ein \glqq Hals\grqq ausbildet, welcher es einem Körper ermöglicht von außerhalb des Mondes ($x>x_2$) zur inneren Region ($|x|<1-\mu-R_2$, wobei $R_1=$Radius der Erde, $R_2=$Radius des Mondes) zu gelangen und umgekehrt. Ein umkreisen der Erde ist aufgrund der Energiebilanz $\Omega_2<\Omega_1$ nie ausgeschlossen. Diskutiert man jedoch die Frage: Existiert ein $E$ für die ein Raumfahrzeug von der Erde aus startend den Mond umkreisen kann, ohne das System Erde-Mond verlassen zu können?, so ist dies durchaus der Fall. Die Bedingung dafür lautet $\Omega_1<E<\Omega_2$. \\
Trotz dieser Betrachtungen ist anzumerken, dass energetische Bedingungen nicht genügen, um die Trajektorie des Raumfahrzeugs oder lediglich Charakteristika, wie z.B. Anzahl der Umdrehungen des zu betrachtenden Körpers bzw. das Verhalten des Körpers für $t\to \infty$ vorherzusagen. Die Dynamik des Systems ist durch die Anfangsbedingungen determiniert, wobei sie bei diesem Problem sensitiv von ihnen beeinflusst wird.
\section{Iterationsverfahren}
\subsection{Runge-Kutta-Verfahren}
In Kapitel \ref{ch:grd} ergab sich die Differenzialgleichung \eqref{dfglg}. Aus Gleichung \eqref{glgbed} ist $\vec{f}(z_1,z_2,z_3,z_4)$ bekannt. Folglich betrachten wir das Problem:
\begin{align}\label{grd_glg}
\dot{\vec  z}=\vec{f}(\vec{z})
\end{align}
Aufgrund der niedrigen Konvergenzordnung $q=1$ verwenden wir das Iterationsverfahren \eqref{glgbed} nicht. Stattdessen nutzen wir das folgende Verfahren mit Konvergenzordnung 3:
\begin{align}\label{grd_glg_it}
&\vec{z}^{\;n+1}= \vec{z}^{\;n} + h \left( \frac{1}{6} \vec{k}_1 +  \frac{4}{6}\vec{k}_2 + \frac{1}{6} \vec{k}_3 \right) \\
&\text{mit}\notag \\
&\left\{
\begin{array}{l l l}
\vec{k}_1=\vec{f}\left(\vec{z}^{\;n}\right)\\
\vec{k}_2=\vec{f}\left(\vec{z}^{\;n}+\frac{h}{2} \vec{k}_{1}\right)\\
\vec{k}_3=\vec{f}\left(\vec{z}^{\;n}-h \vec{k}_{1}+2h \vec{k}_{2}\right)
\end{array}
\right.\ \notag
\end{align}
\subsection{Symplektisches Integrationsverfahren}
In der Hamiltonischen Mechanik bietet es sich an das Symplektische Integrationsverfahren anzuwenden. Dabei wird der Gesamthamiltonoperator in zwei Operator aufgespalten. Speziell für unser Problem ergibt sich:

\begin{align}
H = H_1 + H_2 \quad H_1 = \dfrac{(p_x + y)^2}{2}+\dfrac{(p_y-x)^2}{2} \quad H_2 = -\Omega
\end{align}

Aus diesen beiden Operatoren $H_1, H_2$ bekommen wir die Evolutionsoperatoren unter Anwendung der Vorschrift zur Berechnung der Bewegungsgleichungen.

\begin{align*}
\Psi _1 =
&\left(
\begin{array}{c}
\vec{x} + h \nabla _{\vec{p}} H_1 \\
\vec{p}
\end{array}
\right) \\
\Psi _2 =
&\left(
\begin{array}{c}
\vec{x} \\
\vec{p} + h \nabla _{\vec{x}} \Omega
\end{array}
\right)
\end{align*}

Hiermit gelangen wir zu der vollständigen Iterationsvorschrift für dieses Verfahren, dabei werden Halbschritte bei der Ausführung von $ \Psi _2 $ gemacht.
\begin{align*}
\Psi _2 \left(\frac{h}{2}\right) \Psi _1 (h) \Psi _2 \left(\frac{h}{2}\right) \Leftrightarrow
\begin{array}{l}
\vec{p}^{\;n+\frac{1}{2}} = \vec{p}^{\;n} + \frac{h}{2}  \nabla _{\vec{x}}\: \Omega (\vec{x}^{\;n})\\
\vec{x}^{\;n+1} = \vec{x}^{\;n} + h \nabla _{\vec{p}} \: H_1 (\vec{p}^{\;n},\vec{x}^{\;n})\\
\vec{p}^{\;n+1} = \vec{p}^{\;n+\frac{1}{2}} + \frac{h}{2}  \nabla _{\vec{x}}\: \Omega (\vec{x}^{\;n})
\end{array}
\end{align*}



\subsubsection{Energieerhaltung}
Für den Test der Energieerhaltung haben wir beide numerische Verfahren benutzt, um ein Vergleich zu haben. Zur Berechnung der Energie im Symplektischen Verfahren müssen die berechneten Impulse wieder in die entsprechenden Geschwindigkeiten transformiert werden. Die Transformationsgleichungen lauten:

\begin{align}\label{eq:p-v-trafo}
\left\{
\begin{array}{l l }
	p_x &= \dot{x}-y \\
	p_y &= \dot{y}+x
\end{array}
\right.\
\end{align}
\newpage
\begin{figure}[htp!]
	\input{star_energy.tex}
		\centering
		\caption{Energie-Zeit Diagramm zur Überprüfung der Energieerhaltung mit den Werten \\ $v=2.55541$ und $\alpha=1.44513$.}% \mathtt{rad}
	\label{fig:star-energy}
\end{figure}
 Das Symplektische Verfahren zeigt bei allen gewählten Zeitschritten einen Anstieg der Gesamtenergie an. Ein maximalen Anstiege scheinen periodisch aufzutreten. Beim Runge-Kutta-Verfahren hingegen bleibt die Energie erhalten. Eine mögliche Erklärung dafür ist, dass das Runge- Kutta-Verfahren aufgrund der höheren Konvergenzordnung genauer arbeitet. Die Periodizität der großen Energieverluste könnte ein Indiz dafür sein, dass die Trajektorie zu diesen Zeitpunkten die stärkste Richtungsänderung aufweist, wodurch sich somit auch der größte Rechenfehler ergibt. Daraus sieht man wie stark die Wahl des Iterationsverfahren die Flugbahn beeinflussen kann, zumal die Dynamik sensitiv ist.

\section{Interessante Trajektorien}
Allgemein ist das Ende der Flugbahn erreicht, sobald eine der folgenden 4 Bedingungen zutreffen:
\begin{align*}
&\mathrm{if} (|x|>2)\\
&\mathrm{if} (|y|>2)\\
&\mathrm{if}\left((x+\mu)^2+y^2<R_1^2\right) \\
&\mathrm{if}\left(((x-(1-\mu))^2+y^2<R_2^2\right) \\
\text{wobei } R_1&=0.2=\text{Radius der Erde und }R_2=0.01=\text{Radius des Mondes }
\end{align*}
Weiterhin ist aus Gleichung \eqref{eq:p-v-trafo} (S. \pageref{eq:p-v-trafo}) bekannt, wie sich die initialen Impulse  $\vec p=(p_x,p_y)$ aus den initialen Geschwindigkeiten $\vec v=(v_x,v_y)$ berechnen. Die Vorgabe der jeweiligen $\vec v$ erfolgt in den folgenden jeweiligen Abschnitten.
\subsection{Erde-Mond-Erde}\label{unter}
Als erstes untersuchen wir hier Trajektorien denen es, wie im Abschnitt \ref{kap_2} beschrieben, energetisch möglich ist die Erde zu verlassen und den Mond zu umkreisen, jedoch nicht das Erde-Mond-System zu verlassen. Die Flugbahnen enden somit auf dem Mond oder auf der Erde, wobei wir unsere Betrachtung auf die Bahnen beschränken, welche auf der Erdoberfläche enden. Weiterhin ist $j$ die Anzahl der Mondumdrehungen. Diese wird berechnet, indem man annimmt, dass zu Beginn der Flugbahn $j=0$. Weiterhin wird $j$ gemäß der Bedingung \eqref{counter} bis zum Ende der Trajektorie erhöht:
\begin{flalign}\label{counter}
&\mathrm{if}(y(t_n)>0)\{\notag\\
&\quad\mathrm{if}(\mathrm{sgn}(x(t_n)-(1-\mu))\neq\mathrm{sgn}(x(t_{n+1})-(1-\mu)))\quad	j=j+1;\\
&\}\notag\\
&\text{wobei } \mathrm{sgn}(x)\text{ im Abschnitt \ref{mathe} (S. \pageref{mathe}) definiert wurde}\notag
\end{flalign}
Aus der Gleichung \eqref{Energie} ergibt sich der Bereich für die initialen Geschwindigkeitsbeträge $|\vec {v}|$, die diese Trajektorien begünstigen. Ein Korrekturwert von 0.1 wird aufgrund der Dissipativität des Iterationsverfahrens hinzugefügt.
\begin{align}\label{gesch1}
2.60861\approx\sqrt{2\cdot(\Omega_E-\Omega_2)}+0.1>|\vec {v}|\ge\sqrt{2\cdot(\Omega_E-\Omega_1)}\approx2.49541
\end{align}
Hierbei ist $\Omega_E=\Omega(x_E,y_E)$, wobei $(x_E,y_E)$ ein Punkt auf der Erdoberfläche ist. Hierfür wählten wir $(x_E,y_E)=(0.15,0)$. Weiterhin wird die Orientierung der initialen Geschwindigkeit durch den Winkel $\alpha=\arctan{\big(\frac{\dot y}{\dot x}\big)}$ bestimmt, also ist $v_x=|\vec {v}| \cos(\alpha)$ und $v_y=|\vec {v}| \sin(\alpha)$. Somit haben wir zwei Werte, welche wir auf der Suche nach der maximalen Umdrehungszahl um den Mond variieren können $\alpha\text{ und }|\vec {v}|$. Dabei lassen wir $|\vec {v}|$ konstant und suchen den Winkel $\alpha_\text{max}$ mit dem wir die bis dahin meisten Mondumrundungen erzeugt haben. Danach speichern wir den Wert für diesen Winkel und setzen ihn zurück und variieren $|\vec {v}|$. Im Anschluss starten wir wieder den Suchlauf nach $\alpha_\text{max}$. Das Tupel $(|\vec {v}|,\alpha)$ mit den größtem $j$ und dem von uns gewünschten Ende (Weltall oder Erdoberfläche) wird abgespeichert. Die jeweiligen Variationsalgorithmen sind gleich und im Abschnitt \ref{vari} (S. \pageref{vari}) beschrieben. Hierbei sei $m$ die Dimension des Vektors, welcher für den Variationsalgorithmus benötigt wird. Für den Winkel gilt $\alpha\in [0,\frac{\pi}{2})$.\\
\begin{figure}[h]
\input{star_1.tex}
	\centering
		\caption{Simuliert mit dem Symplektischen Integrationsverfahren mit den Werten $m=20$,\\ $\epsilon=10^{-5}$, $|\vec {v}|=2.52937$ und $\alpha=0.229572$.}
\label{fig:star1}
\end{figure}\\
In Abbildung \ref{fig:star1} ist die vierfache Umdrehung des Mondes mit anschließendem Zurückkehren zur Erde dargestellt. Neben der Abbruchbedingung \eqref{Abbruch} (Seite \pageref{Abbruch}), welche hier lediglich für die Variation des Winkels benötigt wurde, musste der Abbruch der Variation beider Variablen $|\vec {v}|$ und $\alpha$ erfolgen sobald $j=4$.
\newpage
\begin{figure}[h]
\input{star_2.tex}
	\centering
		\caption{Simuliert mit dem Symplektischen Integrationsverfahren mit den Werten $m=120$,\\ $\epsilon=10^{-8}$, $|\vec {v}|=2.51215$ und $\alpha=0.177274$.}
\label{fig:star2}
\end{figure}
In Abbildung \ref{fig:star2} ist die achtzehnfache Umdrehung des Mondes mit anschließendem Zurückkehren zur Erde dargestellt. Dies ist die maximale Anzahl an Mondumrundungen mit anschließender Erdrückkehr, welche wir mit dem Symplektischen Integrationsverfahren gefunden haben.\\
\begin{figure}[h]
\input{star_3.tex}
	\centering
		\caption{Simuliert mit dem Runge-Kutta-Verfahren mit den Werten $m=10$, \\$\epsilon=10^{-6}$, $|\vec {v}|=2.4983$ und $\alpha=1.42628$.}
\label{fig:star3}
\end{figure}\\
In Abbildung \ref{fig:star3} ist die 18-fache Umdrehung des Mondes mit anschließendem Zurückkehren zur Erde dargestellt. Dies ist die maximale Anzahl an Mondumrundungen mit anschließender Erdrückkehr, welche wir mit dem Runge-Kutta-Verfahren gefunden haben.
\newpage
Vergleicht man die Trajektorien Abb. \ref{fig:star3} und Abb.\ref{fig:star2}, so erkennt man, dass die mit dem Runge-Kutta-Verfahren simulierten Trajektorien erst eine Erdumdrehung verrichten bevor sie den Mond umkreisen, was bei der Symplektischen Integration nicht der Fall ist. Außer, dass also sich die $\alpha$ stark unterscheiden, benötigt das Runge-Kutta-Verfahren für maximale Mondumdrehungen kleinere Anfangsgeschwindigkeiten als das Symplektische Integrationsverfahren. Unter Betrachtung der Abbildungen \ref{fig:int1} und \ref{fig:int2} (Abschnitt \ref{interessant}, S. \pageref{interessant}) wird ein weiterer Unterschied ersichtlich. Untersucht wurden die Trajektorien die auf der Mondoberfläche enden mit $|\vec {v}|$, wie bereits definiert, als Anfangsgeschwindigkeit. In Abb. \ref{fig:int1} ist mit 17 die maximale Anzahl von Mondumdrehungen mithilfe des Symplektische Integrationsverfahren gefunden worden. In Abb. \ref{fig:int2} wurde mit 1378 eine vergleichsweise große Anzahl an Mondumdrehungen mit dem Runge-Kutta-Verfahren simuliert, welche nicht die maximale Anzahl ist. Das Problem hierbei war die Rechenzeit, welche für eine so große Anzahl an Umdrehungen sehr groß wird.
\subsection{Erde-Mond-Weltall}
Die nächsten für uns interessanten Trajektorien sind diejenigen, welche nach ihren Mondumrundungen das Erde-Mond-System verlassen. Für die initiale Geschwindigkeit $\vec {v}$ gilt:
\begin{align}\label{gesch3}
0.697826=\sqrt{2\cdot(\Omega (x_2+0.1,0)-\Omega_3)}+0.1>|\vec {v}|\ge\sqrt{2\cdot(\Omega_E-\Omega_2)}=0.250861
\end{align}
Die Bedingung \eqref{gesch3} ermöglicht den Trajektorien das Erde-Mond-System um den Lagrangepunkt $L_2$ zu verlassen, wobei es widerum unmöglich ist das System um $L_3$ bzw. an anderen Orten $(x,y)$ zu verlassen. Trotzdessen ist ein Ende der Trajektorie auf der Erd- bzw. Mondoberfläche nicht ausgeschlossen.
\\
\begin{figure}[h]
\input{star_v1.tex}
	\centering
		\caption{Simuliert mit dem Symplektischen Integrationsverfahren mit den Werten $m=20$, \\$\epsilon=10^{-5}$, $|\vec {v}|=2.5246$ und $\alpha=0.275675$.}
\label{fig:starv1}
\end{figure}\\
In Abbildung \ref{fig:starv1} ist die vierfache Umdrehung des Mondes mit anschließendem Zurückkehren zur Erde dargestellt. Die Abbruchbedingung wurde auch hier, wie im Unterabschnitt \ref{unter} beschrieben, modifiziert.
\newpage
\begin{figure}[h]
\input{star_v2.tex}
	\centering
		\caption{Simuliert mit dem Symplektischen Integrationsverfahren mit den Werten $m=120$, \\$\epsilon=10^{-8}$, $|\vec {v}|=2.5206$ und $\alpha=0.119901$.}
\label{fig:starv2}
\end{figure}
In Abbildung \ref{fig:starv2} ist die neunfache Umdrehung des Mondes mit anschließendem Verlassen des Erde-Mond-Systems dargestellt. Dies ist die maximale Anzahl an Mondumrundungen ohne zwischenzeitlicher Rückkehr zum Erdorbit, welche wir mit dem Symplektischen Integrationsverfahren gefunden haben.\\
\begin{figure}[h]
\input{star_v3.tex}
	\centering
		\caption{Simuliert mit dem Runge-Kutta-Verfahren mit den Werten $m=30$,\\ $\epsilon=10^{-6}$, $|\vec {v}|=2.58559$ und $\alpha=0.578022$.}
\label{fig:starv3}
\end{figure} \\
In Abbildung \ref{fig:starv3} ist die 12-fache Umdrehung des Mondes mit anschließendem Verlassen des Erde-Mond-Systems dargestellt. Dies ist die maximale Anzahl an Mondumrundungen ohne zwischenzeitlicher Rückkehr zum Erdorbit, welche wir mit dem Runge-Kutta-Verfahren gefunden haben.\\
\subsection{Weltall-Mond-Erde-Weltall}
Hierbei soll die im Weltall beginnende Trajektorie mehrere Mondumdrehungen aufweisen, bevor sie nach mehreren Erdumrundungen im Weltall verschwindet. Als Startpunkt wählen wir $(x_2+0.1,0)$, wobei $x_2$ die $x$-Koordinate des Lagrangepunktes $L_2$ ist. Für die Suche passen wir hier den Variationsalgorithmus an, indem wir ein $j_E$ definieren, was die Anzahl der Erdumrundungen angibt. Das für $j_E$ verwendete Zählverfahren gleicht dem für $j$ \eqref{counter}. Damit erweitern wir den Variationsalgorithmus um folgende Abbruchbedingung:\\
\begin{align}\label{Abb2}
&\mathrm{if} (j\ge 2)\{\notag\\
&\quad	\mathrm{if}(j_E\ge 2) \text{ Abbruch};\\
&\}\notag
\end{align}
Der nun durch die Bedingung \eqref{Abb2} erweiterte Variationsalgorithmus findet eine Trajektorie, welche mindestens 2 Mond- bzw. 2 Erdumdrehungen aufweist, bevor die Trajektorie im Weltall endet. Weiterhin gestattet dieses Suchverfahren auch die Vorgabe von mehr als jeweils 2 Umrundungen.Für den Winkel gilt $\alpha\in [\pi,\frac{3 \pi}{2})$. Für die Geschwindigkeit $\vec {v}$ gilt weiterhin Bedingung \eqref{gesch3}:\\
\begin{figure}[h]
\input{star_a3.tex}
	\centering
		\caption{Simuliert mit dem Symplektischen Integrationsverfahren mit den Werten $m=100$,\\ $\epsilon=10^{-6}$, $|\vec {v}|=0.557495$ und $\alpha=4.20377$.}
\label{fig:star_a3}
\end{figure} \\
Die in Abb. \ref{fig:star_a3} dargestellte Trajektorie weist 7 Erdumrundungen und 3 Mondumrundungen auf, wobei eine Mondumrundung erst nach den Erdumkreisungen stattfindet, woraufhin der Körper ins Weltall verschwindet. Bei Vorgabe höherer Umkreisungen sind häufigere Erd-Mond-Übergänge nicht ausgeschlossen, weswegen wir uns bei der Vorgabe der Minimalumrundungen mit jeweils 2 zufrieden gegeben haben. Sofern mehrere Umrundungen ohne Erhöhung der Erd-Mond-Transitionen gefordert wären, müsste man den Suchalgorithmus weiterhin modifizieren. \\
Abb. \ref{fig:star_a33} (S. \pageref{fig:star_a33}) zeigt eine Trajektorie mit mehreren Übergängen und jeweils 4 Umrundungen.
\newpage
\part{Mathematische Grundlagen und Algorithmen}
\section{Bisektionsverfahren} \label{mathe}
Sei $f$ eine stetige Funktion, wobei $f(x_0)$ und $f(y_0)$ verschiedene Vorzeichen haben, so existiert laut des Zwischenwertsatzes eine Nullstelle $z \in [x_0;y_0]$. Durch folgende Iteration kann man $z$ approximativ bestimmen.\\
\begin {align}\label{gleichung bisec}
\left\{
\begin{array}{l l}
z_k=\frac{x_{k-1}-y_{k-1}}{2}, \quad \forall \; k \in \mathbb{Z^+}\\
x_k=z_k,y_{k}=y_{k-1}, \quad \text{falls} \quad \mathop{\mathrm{sgn}}(f(z_k))=\mathop{\mathrm{sgn}}(f(x_{k-1}))\\
y_k=z_k,x_{k}=x_{k-1}, \quad \text{sonst}
\end{array}
\right.\
\end{align}\\
In der Gleichung \eqref{gleichung bisec} ist $\mathop{\mathrm{sgn}}(x)$ die Signumfunktion für die folgendes gilt:
\begin{align*}
\mathop{\mathrm{sgn}}(x)=
\left\{
\begin{array}{l l}
1, \quad \forall \; x>0\\
0, \quad  \text{falls} \quad x=0\\
-1, \quad \text{sonst}
\end{array}
\right.\
\end{align*}\\
\section{Variationsalgorithmus} \label{vari}
Sei $\vec x(t)$ eine Trajektorie um den Mond, wobei die Anzahl der Mondumdrehung $j$ maximiert werden soll, indem der Anfangswert $w$ variiert wird. Weiterhin soll gelten $w \in [a,b)$. Man definiert $\vec W\in \mathbb{R}^m$ und $\vec J\in \mathbb{N}^m$, wobei $\vec W=(w_0, w_1,..,w_{m-1})$ und $\vec J=(j_0, j_1,..,j_{m-1})$. Dabei ist $w_i=a+i\cdot w_{\mathrm{step}} \; \forall \; i\in\{0,1,2,..,m-1\} \text{ mit } w_{\mathrm{step}}=\frac{b-a}{m}$ und $j_i$ gibt die Anzahl an Umdrehungen an, welche mit dem Anfangswert $w_i$ erreicht wurden. Definiert man nun:
\begin{align*}
j_{\mathrm{max}}&=\max\{j_0, j_1,...,j_{m-1}\}\\
l&=\text{Anzahl der Maximalelemente in } \vec J\\
k&=\text{kleinster Index aller Maximalelemente in } \vec J
\end{align*} \\
Das heißt, sofern es $l>1$ Maximalelemente in $\vec J$ gibt, ist $k$ das Minimum der Menge der Indizes aller $j_{\mathrm{max}}$ in $\vec J$.\\
 Mithilfe dieser Grundlagen lautet der Variationsalgorithmus:
\begin{flalign}
&\mathrm{if}(l==m)\notag\\
&\quad w_\mathrm{step, neu}=\epsilon ;\notag\\
&\mathrm{else}\{\notag\\
&	\quad\mathrm{if}(l==1)\{\notag\\
&\quad\quad		\mathrm{if} (k==0)\{\notag\\
&\quad \quad\quad w_0=w_{k};\notag\\
&\quad	\quad\quad		w_\mathrm{step, neu}=\frac{ w_\mathrm{step, alt}}{m};\notag\\
&\quad\quad		\}\notag\\
&\quad\quad		\mathrm{if}(k==m-1)\{\notag\\
&\quad \quad\quad w_0=w_{k-1};\notag\\
&\quad	\quad\quad		w_\mathrm{step, neu}=\frac{ w_\mathrm{step, alt}}{m};\notag\\
&\quad\quad			\}\label{algo}\\
&\quad\quad\mathrm{else}\{\notag\\
&\quad \quad\quad w_0=w_{k-1};\notag\\
&\quad	\quad\quad		w_\mathrm{step, neu}=2\cdot\frac{ w_\mathrm{step, alt}}{m};\notag\\
&\quad\quad\}\notag\\
&	\quad	\}\notag\\
&	\quad	\mathrm{else}\{\notag\\
&	\quad	\quad	w_0=w_k;\notag\\
&\quad	\quad		w_\mathrm{step, neu}=(l-1)\cdot\frac { w_\mathrm{step,alt}}{m};\notag\\
&\quad		\}\notag\\
&	\}\notag
\end{flalign}
Wobei dieser Algorithmus endet, sofern die Abbruchbedingung \eqref{Abbruch} erfüllt ist.
\begin{align}\label{Abbruch}
\mathrm{if}(w_\mathrm{step, neu}\le\epsilon)
\end{align}
Sofern alle in $\vec W$ enthaltenen Anfangswerte die gleiche Anzahl an Mondumrundungen zur Folge haben, erfolgt der Abbruch. Existieren mehrere aber weniger als $m$ Maximalwerte, so beginnt der Algorithmus erneut ab dem Anfangswert $j_\mathrm{max}$ mit dem kleinsten Index $k$ und endet bei $j_\mathrm{max}$ mit dem größten Index. Es wird also nur zwischen diesen $j_\mathrm{max}$ geschaut, wobei wir davon ausgegangen sind, dass diese Maximalwerte hintereinanderliegen, also dass, wenn $l$ die Anzahl der $j_\mathrm{max}$ ist, der größter Index der Maximalwerte $k+l-1$ ist. Folglich ist $w_\mathrm{step, neu}=(l-1)\cdot\frac { w_\mathrm{step,alt}}{m}$. Tritt der Fall auf, dass es genau ein $j_\mathrm{max}$ gibt und $k$ dessen Index sei, so wird um diesen Maximalwert geschaut, also zwischen $w_{k-1}$ und $w_{k+1}$. Folglich ist hier $w_\mathrm{step, neu}=2\cdot\frac{ w_\mathrm{step, alt}}{m}$. Da dies für die Fälle, dass $k=0$ bzw. $k=m-1$ nicht möglich ist, wird hier zwischen $w_{0}$ und $w_{1}$ bzw. $w_{m-2}$ und $w_{m-1}$ geschaut. Somit ist die neue Schrittweite $w_\mathrm{step, neu}=\frac{ w_\mathrm{step, alt}}{m}$.
\part{Anhang}
\section{Andere interessante Trajektorien}\label{interessant}
\begin{figure}[h]
\input{int1.tex}
	\centering
		\caption{Simuliert mit dem Symplektischen Integrationsverfahren mit den Werten $m=120$,\\ $\epsilon=10^{-8}$, $|\vec {v}|=2.54291$ und $\alpha=0.235626$.}
\label{fig:int1}
\end{figure}
\begin{figure}[h]
\input{int2.tex}
	\centering
		\caption{Simuliert mit dem Runge-Kutta-Verfahren mit den Werten $m=10$,\\ $\epsilon=10^{-6}$, $|\vec {v}|=2.4983$ und $\alpha=1.40115$.}
\label{fig:int2}
\end{figure}
\newpage
\begin{figure}[h]
\input{star_a33.tex}
	\centering
		\caption{Simuliert mit dem Symplektischen Integrationsverfahren mit den Werten $m=100$,\\ $\epsilon=10^{-6}$, $|\vec {v}|=0.56685$ und $\alpha=4.21256$.}
\label{fig:star_a33}
\end{figure}
\section{Programmiercode}
%\setlength{\marginparwidth}{2cm}
%\marginpar{\color{red}Entscheiden welches Format wir benutzen wollen}

\begin{comment}
\lstdefinestyle{customc}{
	language=c++,showspaces=false,showstringspaces=false,breaklines=true,numbers=left,frame=single,extendedchars=true,inputencoding=latin1,literate=%
	{Ö}{{\"O}}1
	{Ä}{{\"A}}1
	{Ü}{{\"U}}1
	{ß}{{\ss}}2
	{ü}{{\"u}}1
	{ä}{{\"a}}1
	{ö}{{\"o}}1,
	basicstyle=\footnotesize\ttfamily,
	keywordstyle=\bfseries\color{green!40!black},
	commentstyle=\itshape\color{purple!40!black},
	identifierstyle=\color{blue},
	stringstyle=\color{orange},
	numberstyle=\tiny\color{gray},
}
\lstinputlisting[style=customc]{star_trek.cpp}
\end{comment}

\inputminted[linenos=true,autogobble=true,frame=single,numberblanklines=true,showspaces=false,breaklines=true,breakanywhere=true,framesep=2mm,
baselinestretch=1.2]{cpp}{star_trek.cpp}
\end{document}



```

*make.sh*
```bash
noweb.py -Rtest.tex test.md > test.tex && lualatex -interaction=nonstopmode -shell-escape test.tex && lualatex -interaction=nonstopmode -shell-escape test.tex && notify-send -a Compilation  2025-08-28 fertig && xdg-open test.pdf 2>/dev/null 
```



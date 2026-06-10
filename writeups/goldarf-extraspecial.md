# Extraspecial reframing of the Gold-quadric problem (goldarf.tex draft)

> Moved here from the former `BRIDGES-DRAFT.md`. This is **paste-ready LaTeX** for
> `writeups/goldarf.tex` (preamble additions, a new section, and bibliography
> entries), plus integrator notes. Claim levels are stated inline: the
> extraspecial-group lemmas are standard math with self-contained proofs, the `R_8`
> kernel corollary is implemented and tested (`experiments/misere_kernel.py`),
> reading `E`-equivariance as the naturality criterion is interpretation, and the
> existence of a game-native model of the extension is open (`OPEN.md`).

% ============================================================================
% PASTE 1 — preamble additions for writeups/goldarf.tex
% (place with the existing \newtheorem declarations; the file already has
%  observation/proposition/question in plain style and remark in definition
%  style, so add:)
% ============================================================================

\theoremstyle{plain}
\newtheorem{lemma}{Lemma}
\newtheorem{corollary}{Corollary}
\theoremstyle{definition}
\newtheorem{definition}{Definition}

% ============================================================================
% PASTE 2 — new section. Suggested placement: after the subsection "The
% diagonal framing" (end of "What blocks the game semantics"), before
% "Validation status". Uses existing macros \F, \Arf, \Tr and the existing
% label prop:nogo; new citations Quillen71, Gorenstein80, Winter72 are in
% PASTE 3.
% ============================================================================

\section{The extraspecial reframing}

The diagonal-framing question can be restated in group-theoretic terms, and
the restatement sharpens it. The dictionary is classical -- the
characteristic-$2$ Heisenberg/extraspecial picture of quadratic
forms~\cite{Quillen71,Gorenstein80} -- and this section proves the pieces we
use, to keep the note self-contained. Claim levels: Lemmas~\ref{lem:extdict},
\ref{lem:arftypes}, \ref{lem:abelian} and
Proposition~\ref{prop:between} are standard mathematics; the $R_8$
instantiation in Corollary~\ref{cor:kernel} is implemented and tested
(\path{experiments/misere_kernel.py}); reading $E$-equivariance as the right
naturality criterion is interpretation; the existence of a game-native model
of the extension is open.

Throughout, $V=\F_2^n$ written additively, and $Z=\{1,z\}\cong\mathbb{Z}/2$
written multiplicatively. A \emph{quadratic map} $Q\colon V\to\F_2$ is a
function with $Q(0)=0$ whose polarization
$B(x,y)=Q(x+y)+Q(x)+Q(y)$ is bilinear; $B$ is then alternating (hence
symmetric), and $Q$ may have any diagonal $q_i=Q(e_i)$, exactly the
char-$2$ datum the engine keeps separate from $b$.

\subsection{Quadratic forms are central extensions}

\begin{lemma}[extraspecial dictionary]\label{lem:extdict}
Let $1\to Z\to E\xrightarrow{\ \pi\ }V\to 1$ be a central extension of the
elementary abelian $2$-group $V$ by $Z\cong\mathbb{Z}/2$.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item For $x\in V$ and any lift $\tilde x\in\pi^{-1}(x)$ one has
$\tilde x^2\in Z$, and $\tilde x^2$ is independent of the lift; the
\emph{squaring form} $Q_E\colon V\to\F_2$ defined by
$\tilde x^2=z^{Q_E(x)}$ is well defined, as is the \emph{commutator pairing}
$B_E$ defined by $[\tilde x,\tilde y]=z^{B_E(x,y)}$, which is bilinear and
alternating. They satisfy
\[
(\tilde x\tilde y)^2=\tilde x^2\,\tilde y^2\,[\tilde x,\tilde y],
\qquad\text{equivalently}\qquad
Q_E(x+y)=Q_E(x)+Q_E(y)+B_E(x,y),
\]
so $Q_E$ is a quadratic map with polar form $B_E$.
\item Conversely, every quadratic map $Q$ arises: fixing a basis
$e_1,\dots,e_n$ of $V$, the bilinear $2$-cocycle
\[
\varphi(x,y)=\sum_i x_iy_i\,Q(e_i)+\sum_{i>j}x_iy_j\,B(e_i,e_j)
\]
defines a group $E_Q=Z\times V$ with
$(s,x)(t,y)=(s+t+\varphi(x,y),\,x+y)$, whose squaring form is $Q$ and whose
commutator pairing is the polar form $B$.
\item The assignment $Q\mapsto[E_Q]$ is a bijection from quadratic maps on
$V$ to equivalence classes of central extensions of $V$ by $\mathbb{Z}/2$.
\item $E_Q$ is abelian iff $B=0$ (iff $Q$ is linear); and for $n\geq 1$,
$E_Q$ is \emph{extraspecial} (i.e.\ $Z(E)=[E,E]=\Phi(E)\cong\mathbb{Z}/2$)
iff $B$ is nondegenerate, which forces $n=2r$ even and $|E_Q|=2^{1+2r}$.
\end{enumerate}
\end{lemma}

\begin{proof}
(i) $\pi(\tilde x^2)=2x=0$, so $\tilde x^2\in Z$; replacing $\tilde x$ by
$z\tilde x$ multiplies the square by $z^2=1$. Commutators of lifts lie in
$Z$ because $V$ is abelian, are central, and are unchanged by central
retagging of the lifts; centrality of commutators gives bimultiplicativity,
and $[\tilde x,\tilde x]=1$ gives $B_E(x,x)=0$. For the identity: with
$c=[\tilde y,\tilde x]$ central,
$\tilde x\tilde y\tilde x\tilde y
=\tilde x(\tilde x\tilde y\,c)\tilde y=\tilde x^2\tilde y^2c$,
and $c=c^{-1}=[\tilde x,\tilde y]$ since $Z$ has exponent $2$. As
$\tilde x\tilde y$ lifts $x+y$, the additive identity follows.

(ii) A bilinear map is a $2$-cocycle, so $E_Q$ is a group with center
containing $Z=\{(0,0),(1,0)\}$. Squares: $(s,x)^2=(\varphi(x,x),0)$ and
$\varphi(x,x)=\sum_i x_iQ(e_i)+\sum_{i>j}x_ix_jB(e_i,e_j)=Q(x)$ by the
polarization expansion of $Q$ in the basis. Commutators:
$(s,x)(t,y)(s,x)^{-1}(t,y)^{-1}=(\varphi(x,y)+\varphi(y,x),0)
=(B(x,y),0)$, using that the symmetrization of $\varphi$ is $B$ (the
diagonal terms cancel and $B$ is symmetric).

(iii) Equivalent extensions have cohomologous cocycles, and a coboundary
$\delta\psi(x,y)=\psi(x)+\psi(y)+\psi(x+y)$ has zero diagonal
($\delta\psi(x,x)=\psi(0)=0$), so the squaring form is an invariant of the
class; (ii) gives surjectivity. Injectivity: if two extensions have equal
squaring form $Q$, choose sections and cocycles $\varphi,\varphi'$; both
have diagonal $Q$, and both have symmetrization equal to the polar form of
$Q$, so $d=\varphi+\varphi'$ is a symmetric cocycle with zero diagonal. The
extension $E_d$ is then abelian of exponent $2$, i.e.\ an $\F_2$-vector
space, so $1\to Z\to E_d\to V\to 1$ splits and $d$ is a coboundary; the two
extensions are equivalent.

(iv) Commutators generate $z^{\,\mathrm{im}\,B}$, so $E_Q$ is abelian iff
$B=0$. In general $Z(E_Q)=\pi^{-1}(\operatorname{rad}B)$, since $\tilde x$
is central iff $B(x,\cdot)=0$; thus $Z(E_Q)=Z$ iff $B$ is nondegenerate. If
so, $B\neq0$ gives $[E,E]=Z$, and $\Phi(E)=E^2[E,E]=Z$ because all squares
lie in $Z$ and $E/Z$ is elementary abelian. A nondegenerate alternating
form exists only in even dimension.
\end{proof}

\begin{remark}
This is the characteristic-$2$ Heisenberg picture: $E_Q$ is the Heisenberg
group of $(V,B)$ and $Q$ selects the central extension; cohomologically,
(iii) is the standard identification of $H^2(V;\F_2)$ with quadratic maps
$V\to\F_2$~\cite{Quillen71}. The same finite-quadratic-module data drives
the Weil $S$/$T$ matrices in the library's integral layer. Note that the
dictionary is exactly the engine's discipline in group clothing: the
diagonal $q_i$ (squares) and the off-diagonal $b_{ij}$ (commutators) are
independent central data, and collapsing them is collapsing $E_Q$ onto an
abelian quotient.
\end{remark}

\subsection{Arf classifies the two extraspecial types}

Write $H$ for the hyperbolic plane ($Q(e_1)=Q(e_2)=0$, $B(e_1,e_2)=1$;
$\Arf=0$) and $A$ for the anisotropic plane ($Q\equiv1$ on $V\setminus\{0\}$;
$\Arf=1$), and $\circ$ for the central product (identify the centers).

\begin{lemma}[Arf = the two types]\label{lem:arftypes}
Let $B$ be nondegenerate, $\dim V=2r\geq2$.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item $E_{Q\perp Q'}\cong E_Q\circ E_{Q'}$.
\item $E_H\cong D_8$ and $E_A\cong Q_8$.
\item Every extraspecial $2$-group of order $2^{1+2r}$ is isomorphic to
$E_Q$ for a nondegenerate $Q$, and
$E_Q\cong E_{Q'}$ iff $(V,Q)\cong(V',Q')$ iff the ranks and Arf invariants
agree~\cite{Dickson01,Arf41}. Hence there are exactly two extraspecial
groups of order $2^{1+2r}$:
\[
E^+_{2^{1+2r}}\;=\;D_8^{\circ r}\ \ (\Arf=0,\ Q\cong H^{\perp r}),
\qquad
E^-_{2^{1+2r}}\;=\;D_8^{\circ(r-1)}\circ Q_8\ \ (\Arf=1,\ Q\cong
H^{\perp(r-1)}\perp A).
\]
In particular $D_8\circ D_8\cong Q_8\circ Q_8$, the group avatar of
$H\perp H\cong A\perp A$ (additivity of $\Arf$).
\item (Census.) With $\varepsilon=(-1)^{\Arf Q}$,
\[
\#\{g\in E_Q: g^2=1\}=2\,\#\{x:Q(x)=0\}=2^{2r}+\varepsilon\,2^{r},
\]
so the element-order census of $E_Q$ is the zero-count bias of
Section~5 read multiplicatively, and it distinguishes the two types.
\end{enumerate}
\end{lemma}

\begin{proof}
(i) The cocycle $\varphi\oplus\varphi'$ on $V\oplus V'$ is bilinear with
diagonal $Q\perp Q'$, and
$(s,(x,x'))\mapsto[((s,x),(0,x'))]$ is an isomorphism onto
$(E_Q\times E_{Q'})/\langle(z,z')\rangle$.

(ii) $E_H$ has order $8$, is nonabelian ($B\neq0$), and has a noncentral
involution (any lift of $e_1$ squares to $z^{Q(e_1)}=1$), so $E_H\cong D_8$;
in $E_A$ every noncentral element squares to $z$ (as $Q\equiv1$ off $0$), so
$z$ is the unique involution and $E_A\cong Q_8$.

(iii) For $E$ extraspecial, $V=E/Z(E)$ is elementary abelian (as
$\Phi(E)=Z(E)$), so Lemma~\ref{lem:extdict} applies to
$1\to Z(E)\to E\to V\to1$ and gives $E\cong E_{Q_E}$ with $B_{Q_E}$
nondegenerate ($Z(E)$ exactly central). If $\psi\colon E_Q\to E_{Q'}$ is an
isomorphism, it preserves the (characteristic) centers, fixes $z\mapsto z'$
(the unique nontrivial central element), and induces an $\F_2$-linear map
$\bar\psi\colon V\to V'$ with
$Q'(\bar\psi x)$ read from $\psi(\tilde x)^2=\psi(\tilde x^2)$ -- an
isometry. Conversely an isometry $g$ transports the cocycle:
$\varphi'\circ(g\times g)$ has diagonal $Q'\circ g=Q$, hence by
Lemma~\ref{lem:extdict}(iii) yields an extension equivalent to $E_Q$, and
$(s,x)\mapsto(s,gx)$ closes the isomorphism. Dickson's classification of
nonsingular quadratic forms over $\F_2$ (two isometry classes per even
rank, separated by $\Arf$) finishes; the central-product normal forms
follow from (i), (ii) and Witt decomposition.

(iv) The two lifts of $x$ both square to $z^{Q(x)}$, so they are
involutions or the identity iff $Q(x)=0$; now apply the zero-count formula
of Section~5. The counts differ for $\varepsilon=\pm1$, so
$E^+\not\cong E^-$.
\end{proof}

\subsection{The abelian obstruction}

\begin{lemma}[abelian obstruction]\label{lem:abelian}
Let $M$ be a commutative monoid, $Z=\{1,z\}\subseteq M$ a two-element
subgroup, $N\subseteq M$ a submonoid containing $Z$, and
$\pi\colon N\twoheadrightarrow V$ a surjective monoid homomorphism onto an
$\F_2$-vector space with $\pi^{-1}(\pi(m))=mZ$ for all $m\in N$. Then the
squaring form $Q\colon V\to\F_2$, defined by $\tilde m^2=z^{Q(x)}$ for any
$\tilde m\in\pi^{-1}(x)$, is well defined and \emph{$\F_2$-linear}; its
polar form vanishes identically. Consequently no commutative monoid
realizes, through its own squaring map, a quadratic form with nonzero --
a fortiori nondegenerate -- characteristic-$2$ polar form. Equivalently, by
Lemma~\ref{lem:extdict}(iv): if $B\neq0$, the extension $E_Q$ is
nonabelian and admits no model $(N,Z,\pi)$ inside any commutative monoid.
\end{lemma}

\begin{proof}
$\pi(m^2)=2\pi(m)=0$, so $m^2\in\pi^{-1}(0)=Z$, and
$(zm)^2=z^2m^2=m^2$, so $Q$ is well defined. Commutativity gives
$(mn)^2=m^2n^2$; since $\tilde m\tilde n$ is a preimage of $x+y$, this is
$Q(x+y)=Q(x)+Q(y)$, i.e.\ $B\equiv0$.
\end{proof}

\begin{corollary}[misère quotients: the kernel shadow is linear]\label{cor:kernel}
Let $\mathcal{Q}$ be a finite misère quotient -- a finite commutative
monoid -- with kernel $K$, the maximal subgroup of $\mathcal{Q}$ (the
mutual-divisibility class of the product of all idempotents)~\cite{PS}.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item Every configuration $(N,Z,\pi)$ as in Lemma~\ref{lem:abelian} inside
$\mathcal{Q}$ has linear squaring form: the intrinsic multiplication of a
misère quotient cannot supply a quadratic refinement with $B\neq0$, on $K$
or anywhere else in $\mathcal{Q}$.
\item Under the regularity hypothesis of \cite[Thm.~6.4]{PS} (satisfied by
the regular finite quotients arising in practice), $K\cong(\mathbb{Z}/2)^k$
is the group of normal-play Grundy values and the $P$-portion meets $K$ in
the XOR-linear normal-play set.
\item On the smallest wild quotient
$R_8=\langle a,b,c\mid a^2=1,\ b^3=b,\ bc=ab,\ c^2=b^2\rangle$ with
$P=\{a,b^2\}$: the idempotents are $\{1,b^2\}$, the kernel is
$K=\{b^2,b,ab,ab^2\}\cong(\mathbb{Z}/2)^2$ with identity $b^2$,
$P\cap K=\{b^2\}\mapsto\{0\}$ is linear, and the genuinely misère
$P$-element $a$ lies outside $K$. (Implemented and tested:
\path{experiments/misere_kernel.py}.)
\end{enumerate}
So the linear obstruction observed on $R_8$ is forced by commutativity, not
an accident of the example: a nondegenerate polar form is the commutator
pairing of a nonabelian group, and a misère quotient has none to offer.
\end{corollary}

\begin{proof}
(i) is Lemma~\ref{lem:abelian} applied to submonoids of $\mathcal{Q}$,
which are commutative. (ii) is quoted from \cite{PS}. (iii) is a finite
verification.
\end{proof}

\subsection{A Tier-2 naturality screen: $E$-equivariance}

Proposition~\ref{prop:nogo} kills rules that are blind to everything but
the symplectic structure (Tier~1), while per-position evaluator circuits
for $Q_a$ realize the quadric tautologically (Tier~3). The extraspecial
dictionary suggests a criterion strictly between the two: the rule should
see the extension $E_Q$ -- which carries the diagonal data $q_i$
structurally, as squares -- but only up to its automorphisms, so that no
basis, field structure, or evaluation circuit can be smuggled in.

\begin{definition}[uniform rules; the three tiers]\label{def:tier2}
Let $Q$ be nondegenerate on $V=\F_2^{2r}$ with polar form $B$, and
$E=E_Q$, $\pi\colon E\to V$, $\Sigma=\{x\in V: Q(x)=0\}$.
\begin{enumerate}[label=(\alph*),leftmargin=*]
\item A \emph{uniform rule on a finite set $X$} is a single move relation
$M\subseteq X\times X$, fixed once for all positions; for normal play we
require the move digraph acyclic, and outcomes are computed by the usual
retrograde recursion. (For loopy or misère readings, replace the outcome
recursion; the equivariance constraint below applies verbatim, since
outcome classes are invariants of digraph automorphisms.)
\item A uniform rule on $E$ is \emph{$E$-equivariant} if every
$\alpha\in\operatorname{Aut}(E)$ is an automorphism of the move digraph.
Its \emph{shadow} is $\pi(P)\subseteq V$, where $P$ is its $P$-set; it
\emph{realizes} $T\subseteq V$ if $\pi(P)=T$.
\item \emph{Tier 1} ($\operatorname{Sp}(B)$-blind): a uniform rule on $V$
invariant under $\operatorname{Sp}(B)$, as in
Proposition~\ref{prop:nogo}. \emph{Tier 3} ($Q_a$-evaluation): a family of
games $\{\Gamma_x\}_{x\in V}$ whose structure varies with the designated
input $x$ (e.g.\ the trace-circuit evaluator); this is not a uniform rule.
\emph{Tier 2 screen}: a route to the Gold quadric is Tier-2 natural
\emph{only if} its positions and moves lift to an $E$-equivariant uniform
rule on $E_{Q_a}$ -- with $Q_a$ restricted to its nonsingular core when the
form is degenerate, matching the classifier's radical discipline -- that
realizes $\{Q_a=0\}$.
\end{enumerate}
\end{definition}

\begin{proposition}[the screen sits strictly between the solved tiers]\label{prop:between}
Let $Q$ be nondegenerate on $V=\F_2^{2r}$, $r\geq2$, and $E=E_Q$.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item The image of $\operatorname{Aut}(E)\to GL(V)$ is exactly the
orthogonal group $O(Q)$, with kernel the central twists
$g\mapsto g\,z^{\ell(\pi g)}$, $\ell\in V^{*}$, which equal
$\operatorname{Inn}(E)\cong V$ because $B$ is nondegenerate; thus
$\operatorname{Out}(E)\cong O(Q)$ (cf.~\cite{Winter72}). The
$\operatorname{Aut}(E)$-orbits on $E$ are
\[
\{1\},\qquad \{z\},\qquad \pi^{-1}(\Sigma\setminus\{0\}),\qquad
\pi^{-1}(V\setminus\Sigma).
\]
\item Hence the $P$-set of any $E$-equivariant uniform rule is a union of
these four orbits, and its shadow is one of the eight unions of $\{0\}$,
$\Sigma\setminus\{0\}$, $V\setminus\Sigma$. The quadric $\Sigma$ is among
them: the Tier-1 exclusion does not apply.
\item (Tier 1 $\subsetneq$ Tier 2.) Every $\operatorname{Sp}(B)$-invariant
uniform rule on $V$ pulls back along $\pi$ to an $E$-equivariant rule with
$P$-set $\pi^{-1}$ of the original, whose shadow is therefore one of
$\varnothing$, $\{0\}$, $V\setminus\{0\}$, $V$ -- never the quadric
(Proposition~\ref{prop:nogo} in this language). The containment is strict:
the \emph{squaring rule}
\[
g\to h\ \text{legal}\iff g^2=z\ \text{and}\ h^2=1
\]
(``move from any order-$4$ position to any position of order at most
$2$'') is $E$-equivariant and acyclic, and its $P$-set is the involution
locus $\{g:g^2=1\}=\pi^{-1}(\Sigma)$, with shadow exactly $\Sigma$.
\item (Tier 2 $\subsetneq$ Tier 3.) Up to the isometry induced by any group
isomorphism, the family of $E$-equivariant rules and their shadows depends
only on $(r,\Arf Q)$ (Lemma~\ref{lem:arftypes}(iii)). In particular the
screen cannot separate Gold forms of equal core rank and Arf invariant, and
an $E$-equivariant rule has no access to $m$, the exponent $a$, the field
multiplication, the coordinate frame $q_i=Q_a(e_i)$, or any evaluation
circuit. Conversely, arbitrary subsets of $V$ -- all realizable by ad hoc
acyclic lookup games and by Tier-3 evaluators (Section~7.1) -- are not
shadows of $E$-equivariant rules once they leave the eight-set list of
(ii).
\end{enumerate}
\end{proposition}

\begin{proof}
(i) Automorphisms preserve $Z=Z(E)$ and fix $z$
($\operatorname{Aut}(\mathbb{Z}/2)=1$), so the induced map preserves the
squaring form: $Q(\bar\alpha x)$ is read from
$\alpha(\tilde x)^2=\alpha(\tilde x^2)=z^{Q(x)}$; the image lies in
$O(Q)$. Surjectivity: an isometry $g$ transports the standard cocycle as in
the proof of Lemma~\ref{lem:arftypes}(iii), producing an automorphism of
$E_Q$ inducing $g$. The kernel consists of maps $g\mapsto g\,z^{\ell(\pi g)}$
with $\ell$ linear (multiplicativity forces additivity of $\ell$); inner
automorphisms are exactly the twists $\ell=B(\bar g,\cdot)$, and
nondegeneracy makes $\bar g\mapsto B(\bar g,\cdot)$ onto $V^{*}$. Orbits:
$O(Q)$ is transitive on $\Sigma\setminus\{0\}$ and on $V\setminus\Sigma$ by
Witt's extension theorem in its characteristic-$2$ quadratic-space form
\cite{Taylor} (both sets are nonempty for $r\geq2$); the two lifts of any
$x\neq0$ are conjugate under conjugation by $\tilde g$ with $B(\bar g,x)=1$;
and $1$, $z$ are fixed by everything.

(ii) Digraph automorphisms preserve the retrograde outcome recursion, so
$P$ is $\operatorname{Aut}(E)$-invariant; apply (i) and project.

(iii) The pullback $M=\{(g,h):(\pi g,\pi h)\in R\}$ of an
$\operatorname{Sp}(B)$-invariant rule $R$ is preserved by
$\operatorname{Aut}(E)$, since the induced action on $V$ lies in
$O(Q)\subseteq\operatorname{Sp}(B)$; it is acyclic when $R$ is, and an easy
induction along the acyclic rank gives
$P(M)=\pi^{-1}(P(R))$. By the transitivity of
$\operatorname{Sp}(B)$ on $V\setminus\{0\}$ for $r\geq2$, $P(R)$ is a union
of $\{0\}$ and $V\setminus\{0\}$. For the squaring rule: automorphisms
preserve element orders, so the rule is $E$-equivariant; positions with
$g^2=1$ are terminal, hence $P$; positions with $g^2=z$ have a move (to
$1$), hence $N$. Its shadow $\Sigma$ satisfies
$1<|\Sigma|=2^{2r-1}+\varepsilon2^{r-1}<2^{2r}-1$ for $r\geq2$, so it is
none of the four Tier-1 shadows.

(iv) If $\psi\colon E_Q\to E_{Q'}$ is an isomorphism (it exists iff ranks
and Arf agree, Lemma~\ref{lem:arftypes}(iii)), then
$M\mapsto(\psi\times\psi)(M)$ is a bijection between $E$-equivariant rules
carrying $P$-sets to $P$-sets and shadows to shadows through the induced
isometry $\bar\psi$; in particular it carries rules realizing $\{Q=0\}$ to
rules realizing $\{Q'=0\}$. The shadow constraint of (ii) excludes all but
eight subsets, while Tier-3 constructions realize any subset; the
containment is proper.
\end{proof}

\begin{remark}[what the screen does and does not settle]
Equivariance is a screen, not a construction. Over the \emph{abstract}
group $E$ the screen is satisfiable -- Proposition~\ref{prop:between}(iii)
exhibits the squaring rule, which is nothing but the extraspecial squaring
map of Lemma~\ref{lem:extdict} read as a one-move relation. That rule
consults only the multiplication of $E$, so it is a $Q$-evaluator exactly
insofar as $E$ itself is fed in by hand. The reframing therefore relocates
the open problem: instead of \emph{find a play rule that computes $Q_a$},
it asks \emph{exhibit a game-native model of the extension $E_{Q_a}$} --
positions built from game values, multiplication from game constructions --
after which the rule layer is canonical. Lemma~\ref{lem:abelian} is the
sharp constraint on that search: no commutative value-world, in particular
no misère-quotient kernel (Corollary~\ref{cor:kernel}), can host $E$ once
$B\neq0$. The noncommutativity must enter from a structurally available
asymmetry, and the one that normal, misère, and partizan play all possess
-- and that the symmetric polar form $B$ discards -- is the
first-/second-player asymmetry of the move relation. Testing the existing
routes (kernel, loopy Draw-set, misère quotient) against the screen --
full $E$-symmetry versus only $\operatorname{Sp}(B)$ or an abelian
quotient -- is recorded as a progress target in \texttt{OPEN.md}. Claim
levels: the proposition is standard mathematics; treating
$E$-equivariance as \emph{the} naturality criterion is interpretation;
the existence of a game-native $E$ is open. (For $r=1$ the screen degenerates:
on the anisotropic plane $\Sigma=\{0\}$ and $O(Q)=\operatorname{Sp}(B)$,
so the statement is vacuous there; the Gold targets of interest have
$r\geq2$ cores.)
\end{remark}

\begin{question}
Does any game-native construction realize the extraspecial extension
$E_{Q_a}$ -- equivalently, produce the diagonal data $q_i=Q_a(e_i)$ as
squares in a noncommutative game-built structure -- without an evaluation
circuit for $Q_a$? By Lemma~\ref{lem:abelian} the source cannot be any
commutative game-value monoid.
\end{question}

% ============================================================================
% PASTE 3 — bibliography additions (alphabetical slots indicated):
%   after \bibitem{Gold68}:
% ============================================================================

\bibitem{Gorenstein80} D.~Gorenstein, \emph{Finite Groups}, 2nd ed.,
Chelsea Publishing, New York, 1980.

%   after \bibitem{PS}:

\bibitem{Quillen71} D.~Quillen, \emph{The mod~$2$ cohomology rings of
extra-special $2$-groups and the spectra of quadratic forms}, Math.\ Ann.\
\textbf{194} (1971), 197--212.

%   after \bibitem{Taylor}:

\bibitem{Winter72} D.~L.~Winter, \emph{The automorphism group of an
extraspecial $p$-group}, Rocky Mountain J.\ Math.\ \textbf{2} (1972),
159--168.

% ============================================================================
% Notes for the integrator (not for the .tex):
% - Sources verified against repo surfaces: OPEN.md (extraspecial paragraph,
%   tier dichotomy, progress targets), experiments/misere_kernel.py (R8
%   presentation, idempotents {1,b^2}, kernel {b^2,b,ab,ab^2}, P∩K={b^2},
%   a ∉ K, PS Thm 6.4 + regularity caveat), writeups/goldarf.tex
%   (prop:nogo label, \F/\Arf/\Tr macros, existing Taylor/PS/Arf41/Dickson01
%   bibitems, zero-count formula in §5, bench list in §7.1).
% - "Section~5"/"Section~7.1" in the prose refer to "Arf as conditional
%   win-bias" and "Test benches" in the current draft; renumber if sections
%   move.
% - Claim-level discipline: all lemmas/propositions are standard math with
%   self-contained proofs; Corollary (iii) is implemented-and-tested;
%   the criterion-as-naturality reading is flagged interpretation; existence
%   of a game-native E is flagged open. No new mathematical claims beyond
%   standard extraspecial-group theory are introduced.
% ============================================================================

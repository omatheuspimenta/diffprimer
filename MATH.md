# Algorithm Description and Mathematical Formulation

## 1. Introduction

DiffPrimer addresses the "Exclusive Region Identification Problem" in comparative genomics. The objective is to identify genomic subsequences $S$ present in a target taxonomic sets $\mathcal{T}$ that satisfy a uniqueness constraint with respect to a background set $\mathcal{B}$, and subsequently design Polymerase Chain Reaction (PCR) primers $P$ that maximize specificity.

The pipeline integrates three computational paradigms:
1.  **Alignment-Free Sequence Comparison**: Using $k$-mer multisets for region discovery.
2.  **Bit-Parallel Approximate String Matching**: Myers' algorithm for global homology filtering.
3.  **Semiglobal Dynamic Programming**: For precise local primer verification.

---

## 2. Exclusive Region Identification (Alignment-Free)

Let $\Sigma = \{A, C, G, T\}$ be the alphabet of DNA nucleotides. A genome $G$ is a sequence $G \in \Sigma^*$.
Let $w$ denote a $k$-mer (substring of length $k$). The multiset of $k$-mers in $G$ is defined as:
$$ \mathcal{K}_G = \{ (w, c_w) \mid w \in \Sigma^k \} $$
where $c_w$ is the count of occurrences of $w$ in $G$.

### 2.1 Target and Background Sets

Let the Reference Genome be $R$. We define the set of **Target Unique $k$-mers**, $\mathcal{U}_R$, subject to an abundance threshold $\tau$ (experimentally set to $\tau=1$ for single-copy markers):
$$ \mathcal{U}_R = \{ w \in \Sigma^k \mid (w, c_w) \in \mathcal{K}_R \land c_w \le \tau \} $$

Let the set of Background Genomes be $\mathcal{B} = \{B_1, B_2, \dots, B_N\}$. The **Background Universe** $\Omega_B$ is the union of all background $k$-mers:
$$ \Omega_B = \bigcup_{i=1}^N \{ w \mid (w, c_w) \in \mathcal{K}_{B_i} \} $$

### 2.2 Diagnostic Marker Set

The set of candidate diagnostic markers $\mathcal{D}$ is the set difference:
$$ \mathcal{D} = \mathcal{U}_R \setminus \Omega_B $$

**Computation:**
Constructing $\Omega_B$ explicitly is memory-prohibitive. We approximate $\Omega_B$ using a **Bloom Filter** or exact Hash Set optimizations, achieving average-case $O(1)$ lookups. The complexity of constructing $\mathcal{D}$ is $O(|R| + \sum |B_i|)$.

**Contig Assembly:**
Continuous regions are assembled from $\mathcal{D}$. Two $k$-mers $w_i, w_j \in \mathcal{D}$ (at positions $p_i, p_j$ in $R$) are merged into a region $S_{region}$ if $p_j = p_i + 1$. This yields a set of candidate amplicons $\mathcal{S} = \{S_1, S_2, \dots\}$.

---

## 3. Specificity Verification Algorithm

For each candidate region $S \in \mathcal{S}$ and designed primer pair $(P_{fwd}, P_{rev})$, we verify specificity against every background genome $B \in \mathcal{B}$. This is a two-stage process.

### 3.1 Stage I: Global Homology Filter (Myers' Bit-Vector Algorithm)

To efficiently identify if $S$ appears in $B$ with edit distance $d < d_{max}$, we employ the **Myers Bit-Vector Algorithm** (Myers, 1999). This algorithm exploits the bit-parallelism of modern processors to compute a column of the Dynamic Programming (DP) matrix in $O(m/w)$ operations, where $m=|S|$ and $w=64$ (machine word size).

**State Variables:**
Let $D_{i,j}$ be the edit distance between $S[1..i]$ and $B[1..j]$.
We define delta variables for the vertical differences in the DP matrix:
$$ \Delta V_{i,j} = D_{i,j} - D_{i-1,j} \in \{-1, 0, +1\} $$
These are encoded using two bit-vectors, $VP$ (Vertical Positive) and $VN$ (Vertical Negative), such that the $i$-th bit is set if $\Delta V_{i,j} = +1$ or $-1$ respectively.

**Recurrence Relations:**
For each character $B[j]$ in the background:
1.  **Pattern Mask**: $X_j = \text{PM}[B[j]]$ (Bitmask of occurrences of char $B[j]$ in $S$).
2.  **Horizontal/Vertical Interaction**:
    $$ D0 = (((X_j \& VP) + VP) \oplus VP) \lor X_j \lor VN $$
    $$ HP = (VN \lor \neg(D0 \lor VP)) $$
    $$ HN = VP \& D0 $$
3.  **State Update**:
    $$ VP' = (HN \ll 1) \lor \neg(D0 \lor HN \lor X_j) \text{ (approximation)} $$
    $$ VN' = (HP \ll 1) $$
    *(Note: Implementation uses full `long::Myers` chaining for patterns where $m > 64$, maintaining state across `u64` blocks).*

**Complexity:** $O(\lceil \frac{m}{64} \rceil \cdot n)$, where $n=|B|$. This provides a $64\times$ speedup over standard Levenshtein calculation.

### 3.2 Stage II: Local Verification (Semiglobal Alignment)

When Stage I identifies a "Hit" (substring $B'$ in $B$ such that $\text{Lev}(S, B') \le d_{max}$), we must verify if the **primers** bind to $B'$. Due to insertions/deletions (indels), the primers' positions in $B'$ are shifted relative to $S$.

We solve this using **Semiglobal Alignment** (Needleman-Wunsch variant with free end-gaps).
Let $S$ be the query (Region) and $W$ be the target window extracted from $B$ around the hit.

**Score Optimization:**
Maximized Score matrix $H$:
$$ H_{i,j} = \max \begin{cases} H_{i-1, j-1} + \sigma(S[i], W[j]) \\ H_{i-1, j} + \gamma_{\text{gap}} \\ H_{i, j-1} + \gamma_{\text{gap}} \end{cases} $$
where $\sigma(a, b) = 1$ if $a=b$ else $-1$, and $\gamma_{\text{gap}} = -1$.
Initialization: $H_{0,j} = 0, H_{i,0} = i \cdot \gamma$ (Semiglobal in target, global in query).

**Coordinate Mapping:**
An optimal alignment trace $\mathcal{T}$ is computed. We map the primer interval $[x_{start}, x_{end}]$ from $S$ to $[y_{start}, y_{end}]$ in $W$ via the trace function $\tau$:
$$ P_{mapped} = \tau(\mathcal{T}, P_{original}) $$
This accounts for all indels between the primer site and the region boundaries.

**Final Specificity Test:**
We compute the local Hamming/Levenshtein distance inside the mapped sites:
$$ \delta_{local} = \text{BoundedLev}(P_{seq}, W[y_{start}..y_{end}]) $$

If $\delta_{local} \le T_{mismatch}$ (default 2 mismatches) for **both** forward and reverse primers, the amplicon is flagged as **NonSpecific** (potential false positive).

---

## 4. Complexity Analysis

1.  **Region Discovery**: $O(N \cdot L)$, linear with respect to total input size.
2.  **Primer Design**: $O(|\mathcal{S}| \cdot m^2)$ via Primer3 (Thermodynamic alignment).
3.  **Specificity Check**:
    -   Worst case: $O(|\mathcal{S}| \cdot \lceil \frac{m}{64} \rceil \cdot L_{total})$.
    -   Average case is heavily optimized by bit-vector filtering, where expensive Stage II alignment is only performed on rare high-similarity hits ($<1\%$ of genome).

## 5. Statistical Scoring

Homology is reported as:
$$ \text{Similarity Score} = \left( 1 - \frac{d_{min}}{L_{amplicon}} \right) \times 100 $$
where $d_{min}$ is the minimum Levenshtein distance found across all background genomes.

This methodology ensures that DiffPrimer minimizes both False Positives (non-specific binding) and simplifies the computational burden of scanning gigobase-sized genomes.

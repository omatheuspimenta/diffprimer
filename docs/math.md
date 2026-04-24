# Algorithm Description and Mathematical Formulation

## Introduction

DiffPrimer addresses the "Exclusive Region Identification Problem" in comparative genomics. The objective is to identify genomic subsequences $S$ present in a target taxonomic set $\mathcal{T}$ that satisfy a uniqueness constraint with respect to a background set $\mathcal{B}$, and subsequently design Polymerase Chain Reaction (PCR) primers $P$ that maximize target specificity.

The pipeline integrates three computational paradigms:  
1. **Alignment-Free Sequence Comparison**: Utilizing $k$-mer multisets for rapid region discovery.  
2. **Bit-Parallel Approximate String Matching**: Employing Myers' algorithm for global homology filtering.  
3. **Semiglobal Dynamic Programming**: Applying sequence alignment for precise local primer verification.  

---

## Exclusive Region Identification (Alignment-Free)

Let $\Sigma = \{A, C, G, T\}$ be the alphabet of DNA nucleotides. A genome $G$ is defined as a sequence $G \in \Sigma^*$.
Let $w$ denote a $k$-mer (a substring of length $k$). The multiset of $k$-mers in $G$ is given by:

$$
\mathcal{K}_G = \{ (w, c_w) \mid w \in \Sigma^k \}
$$

where $c_w$ represents the number of occurrences of $w$ in $G$.

### Target and Background Sets

Let the Reference Genome be $R$. We define the set of **Target Unique $k$-mers**, $\mathcal{U}_R$, subject to an abundance threshold $\tau$ (empirically set to $\tau=1$ for single-copy markers):

$$
\mathcal{U}_R = \{ w \in \Sigma^k \mid (w, c_w) \in \mathcal{K}_R \land c_w \le \tau \}
$$

Let the set of Background Genomes be $\mathcal{B} = \{B_1, B_2, \dots, B_N\}$. The **Background Universe** $\Omega_B$ is the union of all $k$-mers present in the background:

$$
\Omega_B = \bigcup_{i=1}^N \{ w \mid (w, c_w) \in \mathcal{K}_{B_i} \}
$$

### Diagnostic Marker Set

The set of candidate diagnostic markers $\mathcal{D}$ is obtained via set difference:

$$
\mathcal{D} = \mathcal{U}_R \setminus \Omega_B
$$

**Contig Assembly:**
Continuous regions are assembled from the diagnostic marker set $\mathcal{D}$. Two $k$-mers $w_i, w_j \in \mathcal{D}$ (located at positions $p_i, p_j$ in $R$) are merged into a continuous region if they are adjacent ($p_j = p_i + 1$). This process yields a set of candidate amplicons $\mathcal{S} = \{S_1, S_2, \dots\}$.

---

## Specificity Verification Algorithm

For each candidate region $S \in \mathcal{S}$ and its corresponding designed primer pair $(P_{fwd}, P_{rev})$, the pipeline verifies specificity against every background genome $B \in \mathcal{B}$. This verification is structured as a two-stage process.

### Stage I: Global Homology Filter (Myers' Bit-Vector Algorithm)

To efficiently determine whether $S$ appears in $B$ within an edit distance $d < d_{max}$, DiffPrimer employs the **Myers Bit-Vector Algorithm** (Myers, 1999). This approach leverages the bit-parallelism of modern processor architectures to compute a column of the Dynamic Programming (DP) matrix in $O(m/w)$ operations, where $m=|S|$ and $w=64$ represents the machine word size.

**State Variables:**
Let $D_{i,j}$ denote the edit distance between the prefix $S[1..i]$ and the prefix $B[1..j]$.
We define delta variables to represent the vertical differences in the DP matrix:

$$
\Delta V_{i,j} = D_{i,j} - D_{i-1,j} \in \{-1, 0, +1\}
$$

These variables are encoded using two bit-vectors, $VP$ (Vertical Positive) and $VN$ (Vertical Negative). Specifically, the $i$-th bit is set if $\Delta V_{i,j} = +1$ or $-1$, respectively.

**Complexity:** The time complexity is $O(\lceil \frac{m}{64} \rceil \cdot n)$, where $n=|B|$. This yields a roughly $64\times$ speedup over the standard Levenshtein distance calculation.

### Stage II: Local Verification (Semiglobal Alignment)

When Stage I identifies a "Hit" (a substring $B'$ in $B$ such that $\text{Lev}(S, B') \le d_{max}$), it is necessary to verify whether the designed **primers** can hybridize to $B'$. Due to the potential for insertions and deletions (indels), the primers' coordinates within $B'$ may be shifted relative to their original positions in $S$.

This coordinate mapping is resolved using **Semiglobal Alignment** (a Needleman-Wunsch variant allowing free end-gaps). Let the candidate region $S$ serve as the query and $W$ be the target window extracted from $B$ around the homologous hit.

**Score Optimization:**
The alignment maximizes the score matrix $H$:

$$
H_{i,j} = \max \begin{cases} H_{i-1, j-1} + \sigma(S[i], W[j]) \\ H_{i-1, j} + \gamma_{\text{gap}} \\ H_{i, j-1} + \gamma_{\text{gap}} \end{cases}
$$

where the substitution penalty is $\sigma(a, b) = 1$ if $a=b$ and $-1$ otherwise, and the gap penalty is $\gamma_{\text{gap}} = -1$.
The initialization follows semiglobal constraints (global in the query, local in the target): $H_{0,j} = 0$, and $H_{i,0} = i \cdot \gamma_{\text{gap}}$.

**Coordinate Mapping:**
By tracking the optimal alignment trace $\mathcal{T}$, the primer interval $[x_{start}, x_{end}]$ from the reference sequence $S$ is projected onto the target window $W$ as $[y_{start}, y_{end}]$ via the trace mapping function $\tau$:

$$
P_{mapped} = \tau(\mathcal{T}, P_{original})
$$

This robust mapping explicitly accounts for any indels occurring between the primer binding site and the amplicon boundaries.

**Final Specificity Test (Positional Mismatch Penalty):**
Rather than relying on a naive local edit distance, the algorithm applies a **positional weight** $\omega(x)$ to every mismatch, insertion, or deletion relative to the mapped primer sequence $P_{seq}$ of length $L_P$. This mechanism accurately models the biochemical reality that mismatches at the 3' terminal end severely disrupt polymerase extension, thereby conferring high specificity.

Let $x \in [1, L_P]$ denote the nucleotide position along the primer in the 5' to 3' direction. A configurable **penalty array** $A = [p_{default}, p_{-5}, p_{-4}, p_{-3}, p_{-2}, p_{-1}]$ dictates the mismatch severity. The default weight $p_{default}$ applies to the 5' region, while the subsequent discrete values are assigned to the last five bases at the 3' end (with $p_{-1}$ representing the extreme 3' terminal base). The system defaults to $A = [1.0, 3.0, 3.0, 3.0, 3.0, 3.0]$.

The positional weight function is formulated as:

$$
\omega(x) = \begin{cases} A[x - L_P + 5] & \text{if } x > L_P - 5 \text{ (local 3' region)} \\ A[0] & \text{otherwise} \end{cases}
$$

*(Note: For reverse primers, spatial coordinates are automatically adjusted during alignment. This ensures the 3' weighting accurately reflects hybridization constraints along the forward strand.)*

The positional mismatch score $\delta_{pos}$ is aggregated across the alignment traversal for all editing operations:

$$
\delta_{pos} = \sum_{\text{op} \in \{\text{Subst}, \text{Ins}, \text{Del}\}} \omega(\text{pos}(\text{op}))
$$

**Specificity Classification Tags:**

Based on the global similarity and the local positional mismatch scores ($\delta_{pos}$) evaluated against a configurable threshold $T_{mismatch}$ (default 7.0), each candidate amplicon is classified under one of three definitive tags:

1. **`Specific_LowGlobalSim`**: The global sequence similarity of the candidate region against all off-target genomes falls strictly below the global similarity threshold. The amplicon is globally unique and highly specific.
2. **`NonSpecific_Amplification`**: An off-target region exhibits high homology (similarity $\ge$ threshold) AND the local mismatch score $\delta_{pos}$ remains strictly $< T_{mismatch}$ for **both** the forward and reverse primers. This indicates cross-reactivity and a high probability of false-positive amplification.
3. **`Specific_PositionalMismatches`**: While an off-target region exhibits high homology, the local mismatch score $\delta_{pos}$ is $\ge T_{mismatch}$ for at least one of the primers (frequently driven by critical 3' mismatches). Consequently, amplification is biochemically inhibited, ensuring the primer pair remains highly specific.

---

## Statistical Scoring

Overall sequence homology is standardized and reported as a percentage:

$$
\text{Similarity Score} = \left( 1 - \frac{d_{min}}{L_{amplicon}} \right) \times 100
$$

where $d_{min}$ is the absolute minimum Levenshtein distance identified across the entire background genomic dataset.

This rigorous two-stage methodology ensures that DiffPrimer reliably minimizes false-positive amplifications while substantially reducing the computational overhead required to scan gigabase-sized genomes.

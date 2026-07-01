# Intelligent Reflecting Surface (IRS) vs Decode-and-Forward (DF) Relay Under Stochastic Fading (3GPP SCM)

This repository contains the simulation code, generated figures, and manuscript files for a comparative study of **Intelligent Reflecting Surface (IRS)** and **Decode-and-Forward (DF) relay** performance under realistic stochastic channels using the **3GPP Spatial Channel Model (SCM)**.

The work extends deterministic IRS-vs-relay analysis by incorporating fading randomness, Ricean LOS/NLOS control, pathloss exponent variation, and multipath cluster richness through Monte Carlo evaluation.

---

## 📌 Project Overview
![Graphical Abstract](images/graphical_abstract.png)

Modern wireless systems require both high throughput and robust reliability under realistic channel conditions.  
This project evaluates three schemes in a 3GPP UMi downlink setting:

- **SISO direct transmission**
- **Decode-and-Forward (DF) relay**
- **IRS-assisted transmission** (with coherent phase alignment)

### Main objective
To quantify when and why IRS outperforms DF relay under stochastic fading, and derive practical deployment guidelines based on:

- Ricean **K-factor**
- Pathloss exponent **α<sub>PL</sub>**
- Cluster count **L**
- Minimum element threshold **N<sub>min</sub>**

---

## 🧪 Simulation Scope

### System Setup
- Carrier frequency: **3 GHz**
- Bandwidth: **10 MHz**
- Channel model: **3GPP SCM**
- Scenario: **Urban Micro-cell (UMi)**
- Monte Carlo realizations: **1000 per evaluation point**
- Fixed random seed: **42**

### Parameter Sweeps
1. **K-factor sweep**: \( K \in \{0,10,25\} \) dB  
2. **Pathloss exponent sweep**: \( \alpha_{PL} \in \{0.5,1.0,2.0,3.0,4.5\} \)  
3. **Cluster count sweep**: \( L \in \{1,6,20,50,100\} \)  
4. **Nmin sweep**: minimum IRS elements required to outperform DF relay

---

## 📈 Key Findings (Summary)

- IRS with **N = 200** consistently outperforms DF relay in median SE across tested conditions.
- IRS outage performance is more sensitive to LOS quality (**~1.8×** K-factor sensitivity vs DF at outage level).
- IRS advantage is strongest in guided/low-exponent environments (\( \alpha_{PL} \le 2 \)).
- Increasing cluster richness steepens CDFs (better reliability) without large median shift.
- \( N_{min} \) decreases with K-factor, with diminishing returns beyond ~10 dB.

---

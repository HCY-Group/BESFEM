---
layout: default
title: Governing Equations
nav_order: 3
---

# Governing Equations

## Cathode Concentration

$$\frac{\partial C_c}{\partial t}=\frac{1}{\psi_c}\nabla\cdot\left(\psi_c D_c \nabla C_c\right)-\frac{|\nabla\psi_c|}{\psi_c}r_c$$

---

## Cathode Potential

$$\nabla\cdot\left(\psi_c\kappa_c\nabla\phi_c\right)-|\nabla\psi_c|z_-Fr_c=0$$

---

## Electrolyte Concentration

$$\frac{\partial C_e}{\partial t}=\frac{1}{\psi_e}\nabla\cdot\left(\psi_eD_e\nabla C_e\right)+\frac{|\nabla\psi_c|}{\psi_e}\frac{r_ct_-}{\nu_+}+\frac{|\nabla\psi_a|}{\psi_e}\frac{r_at_-}{\nu_+}$$

---

## Electrolyte Potential

$$\nabla\cdot\left[\psi_e(z_+m_+-z_-m_-)FC_e\nabla\phi_e\right]+|\nabla\psi_c|\frac{r_c}{\nu_+}+|\nabla\psi_a|\frac{r_a}{\nu_+}=\nabla\cdot\left[\psi_e(D_- - D_+)\nabla C_e\right]$$

---

## Anode Concentration

$$\frac{\partial C_a}{\partial t}=\frac{1}{\psi_a}\nabla\cdot\left[\psi_aM_a\nabla\left(\frac{\partial f_G}{\partial C_a}-\varepsilon\nabla^2C_a\right)\right]-\frac{|\nabla\psi_a|}{\psi_a}r_a$$

---

## Anode Potential

$$\nabla\cdot\left(\psi_a\kappa_a\nabla\phi_a\right)-|\nabla\psi_a|z_-Fr_a=0$$

---
TwoDimBoussinesq
================

This gives a very brief introduction to how the Two-Dimensional
Boussinesq equations are derive.

3D Model Equations
------------------

The three-dimensional Boussinesq equations is a simplified model from
the Navier-Stokes equations. It is valid in the limit where the density
variations are small compared to the mean density

.. math:: \rho' \ll \overline{\rho}

This is a very good approximation for ocean dynamics and sometimes
appropriate for the atmosphere. In this limit the density is
approximated everywhere with a constant except in the buoyancy term. A
very general way to write it, as used in the MITgcm is,

.. math::

   \begin{aligned}
   \frac{D\vec u}{Dt} +  2 \vec \Omega \times \vec u 
   + \frac{g \rho}{\overline{\rho}} \hat k + \frac{1}{\overline{\rho}} \vec\nabla p 
   & = \frac{1}{\overline{\rho}} \vec\nabla \cdot \vec \tau, \\
    \vec\nabla \cdot \vec u & = 0, \\
   \frac{D S}{Dt}  & = 0, \\
   \frac{D \theta}{Dt} &= \frac{1}{\overline{\rho} c_{pS}} \vec\nabla \cdot \mathcal{F}_\theta, \\
   \rho & = \rho(\theta, S, z).\end{aligned}

2D Model Equations: Vorticity Approach
--------------------------------------

There are many problems that can be studied in a two-dimensional context
and are still interesting. To study stratified flows we assume that the
fields can vary in one horizontal direction, say :math:`x`, and in the
vertical, :math:`z` We allow for motion in the :math:`y` direction but
no changes in that direction. In the 2D context the problem can be set
up using the pressure method (cite something here) or using the
vorticity method. Since TwoDimboussinesq focuses on the latter we focus
on this approach.

Conservation of mass combined with the fact that no fields vary in the
:math:`y` direction, allow us to define a two-dimensional streamfunction
which then determines the velocity field, :math:`(u,w)`,

.. math::

   u = -\partial_z \psi, 
   \quad \mbox{ and } \quad
   w = \partial_x \psi.

The momentum equations can be written in terms of the Jacobian,

.. math::

   \begin{aligned}
   \partial_t u   & =  J(\psi,  \partial_z \psi) + f v  +\nu \nabla^2 u  -  \frac{1}{{\rho}_0} \partial_x p, \\
   \partial_t v  & =   - J(\psi, v)  -f u  +\nu \nabla^2 v , \\
   \partial_t w   & =  - J(\psi, \partial_x \psi )  + b  +\nu \nabla^2 w  -  \frac{1}{{\rho}_0} \partial_z p.\end{aligned}

To form the :math:`y`-component of vorticity,
:math:`\nabla^2 \psi = - \partial_z u + \partial_x w` we compute the
negative :math:`z` derivative of the :math:`u` eqn and add the :math:`z`
derivative of the :math:`w` eqn,

.. math:: \partial_t \nabla^2 \psi = -J(\psi, \nabla^2 \psi) -  f v_z   + \partial_x b + \nu \nabla^2 \nabla^2 \psi

Therefore the governing equations are,

.. math::

   \begin{aligned}
   \partial_t \nabla^2 \psi &= -J(\psi, \nabla^2 \psi)-  f v_z   +  \partial_x b + \nu \nabla^2 \nabla^2 \psi, \\
   \partial_t v  & =   -J(\psi, v)  + f \partial_z \psi  +\nu \nabla^2 v , \\
   \partial_t b & = -J(\psi, b) + \kappa \nabla^2 b.\end{aligned}

If we assume that there is a background linear stratification denoted
with :math:`\overline{b}(z) = N^2 z` and decompose the buoyancy as

.. math:: b = \overline{b}(z) + b'.

We substitute our decomposition into the equations and get,

.. math::

   \begin{aligned}
   \partial_t \nabla^2 \psi &= -J(\psi, \nabla^2 \psi)-  f v_z   +  \partial_x b + \nu \nabla^2 \nabla^2 \psi, \\
   \partial_t v  & =   -J(\psi, v)  + f \partial_z \psi  +\nu \nabla^2 v , \\
   \partial_t b & = -J(\psi, b) - N^2 \partial_x \psi + \kappa \nabla^2 b.\end{aligned}



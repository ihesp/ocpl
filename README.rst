==================================
ocpl -- a CESM Ocean-Coupler Component
==================================

The "ocean-coupler" (ocpl) is a cesm ocean component that connects to the standard cesm2.1 coupler.  
Ocpl is not a typical cesm ocean model, rather it implements a second-level coupler that
couples/embeds a regional roms model into a global pop model.

The ocpl takes data from the cesm coupler and sends it to roms and pop, 
manages the coupling between roms to pop, merges their surface states, 
and sends this composite ocean surface state back to the cesm coupler as if 
it were simply another cesm-compatible global ocean model.


Also see the CESM web site for documentation and information:

http://www.cesm.ucar.edu

The CESM Quickstart Guide is available at:

http://escomp.github.io/cesm


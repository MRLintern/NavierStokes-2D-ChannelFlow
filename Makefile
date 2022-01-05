FC = gfortran
FCFLAGS = -O3

.PHONY: all clean

all: NavierStokes-2D-ChanelFlow

NavierStokes-2D-ChanelFlow: NS2D.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	$(RM) NavierStokes-2D-ChanelFlow
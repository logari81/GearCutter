# GearCutter 0.2
# Generation of involute gear tooth profiles
# Copyright (C) 2024-2024 Konstantinos Poulios
#
# This file is part of GearCutter.
#
# GearCutter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import gear_cutter_rs as gc

def test_generate_tooth_profile():

  tool = gc.CuttingToolRack(m_n=16,
                            h_aP0=1.25,
                            alfa_P0=20,
                            h_fP0=1,
                            rho_aP0=0.2)
  X,Y,ID = gc.generate_tooth_profile(tool,
                                     z=19,
                                     beta=0.,
                                     x=0.,
                                     d_a=348.,
                                     n_sim=1000,
                                     n_profile=1200,
                                     n_tip=7)

  np.savetxt('tooth_profile.dat', np.column_stack((X,Y,ID)), fmt='%10.5f %10.5f %3d')

if __name__ == '__main__':
  test_generate_tooth_profile()

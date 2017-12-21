function is_same = compare_single_well_methods()
% Compare the bound energies found by the method in tutorial 2 to the bound
% energies found using the propagator method. Return true if all the bound
% energies found are within 1e-8 of each other.
%
% Author: Jack Hong, 30935134
% Last modified: Sept. 29, 2016

E_bound = single_well;

iter = 1;
for E = tutorial2(0.3,1,10)
    if abs(E - E_bound(iter)) > 1e-8
        is_same = false;
        return
    end
    iter = iter + 1;
end

is_same = true;

end
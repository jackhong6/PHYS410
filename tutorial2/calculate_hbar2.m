function hbar2 = calculate_hbar2()
% Return hbar2 in units of m_e * eV *nm^2
m_e = 9.109384e-31;    % kg
eV = 1.602176e-19;     % J
h_bar = 1.054572e-34;  % J*s

hbar2 = h_bar^2 * (1/m_e) * (1/eV) * (1e18);
end
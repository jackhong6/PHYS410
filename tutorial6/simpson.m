function [ Is ] = simpson( fa, fb, fmid, a, b )
% Return integral using Simpson's rule.
Is = (b-a)/6 * (fa + 4*fmid + fb);

end


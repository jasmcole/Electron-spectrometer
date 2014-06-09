function PointingCalc(d,d2,mmerror, lookup)

load(lookup)

C11 = contour(angles*1e3, Energies, screenpos, d - mmerror);
hold on
C12 = contour(angles*1e3, Energies, screenpos, d + mmerror);
C21 = contour(angles*1e3, Energies, screenpos2, d2 - mmerror);
C22 = contour(angles*1e3, Energies, screenpos2, d2 + mmerror);
hold off

% C1 = contour(angles*1e3, Energies, screenpos, d);
% hold on
% C2 = contour(angles*1e3, Energies, screenpos2, d2);

caxis([1e10 1e11])

C11(:,1) = [];
C21(:,1) = [];
C12(:,1) = [];
C22(:,1) = [];

intersection(:,1) = InterX(C11, C21);
intersection(:,3) = InterX(C12, C22);
intersection(:,2) = InterX(C11, C22);
intersection(:,4) = InterX(C12, C21);

% intersection = InterX(C1, C2);

fill(intersection(1,:), intersection(2,:), 1)
% scatter(intersection(1,:), intersection(2,:))
hold off
xlabel('Beam Pointing / mrad');
ylabel('Beam Energy / MeV');

end




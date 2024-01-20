% clear workspace
clear;

% import data
balloon_df = readtable("Sample High Altitude Balloon Data.txt");

% split on the max altitude
[max_alt, max_alt_index] = max(balloon_df{:, "Alt_m"});
ascent_df = balloon_df(1 : max_alt_index, :);
descent_df = balloon_df(max_alt_index+1 : end, :);

% plot alt vs tmp for the whole flight
% plt = plot(balloon_df{:, "TempC"}, balloon_df{:, "Alt_m"}, ".k", MarkerSize=12);
plt = scatter(balloon_df{:, "Alt_m"}./1000, balloon_df{:, "TempC"}, 12, "k", "filled");
set(gca, FontSize=16, FontName="Times New Roman")
ylabel("Temperature [Celsius]")
xlabel("Altitude [kilometers]")
xlim([0, max_alt/1000])
title("Sample flight of a high altitude balloon")
exportgraphics(gcf, "balloon_flight.png")


% export the ascent alt and ascent temp code
ascent_alt = ascent_df{:, "Alt_m"};
ascent_tempC = ascent_df{:, "TempC"};

fileID = fopen("ascent_data.txt", "w+");
fprintf(fileID, '%13s %10s\n', 'ascent_alt', 'ascent_tempC');
fprintf(fileID, '%10.1f %10.1f\n', ascent_alt, ascent_tempC);
fclose(fileID);

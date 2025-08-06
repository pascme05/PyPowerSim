function num = str2number(inputString)

% Split the string by the comma
parts = strsplit(inputString, ',');

% Convert the parts to numeric values
realPart = str2double(parts{1});
imaginaryPart = str2double(parts{2});

% Create the complex number
num = complex(realPart, imaginaryPart);

end


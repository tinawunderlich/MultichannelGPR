function compileMEX(progPath)
% compileMEX  Compile semiCircle.cc if the MEX binary is missing or outdated.
% by Morten Harms 2026
%
% The OpenMP compiler flag differs by platform:
%   Linux / macOS   -fopenmp  (GCC / Clang with libomp)
%   Windows         /openmp   (MSVC)

    srcFile = fullfile(progPath, 'semiCircle.cc');
    mexFile = fullfile(progPath, ['semiCircle.' mexext]);

    needsCompile = ~exist(mexFile, 'file');

    if ~needsCompile
        srcInfo = dir(srcFile);
        mexInfo = dir(mexFile);
        needsCompile = srcInfo.datenum > mexInfo.datenum;
    end

    if needsCompile
        %fprintf('Compiling semiCircle.cc... ');
        if ispc
            % Windows (MSVC)
            mex('-outdir', progPath, ...
                'COMPFLAGS=$COMPFLAGS /openmp', ...
                srcFile);
        else
            % Linux and macOS (GCC / Clang)
            try
                mex('-outdir', progPath, ...
                    'CXXFLAGS=$CXXFLAGS -fopenmp', ...
                    'LDFLAGS=$LDFLAGS -fopenmp', ...
                    srcFile);
            catch
                % Mac:
                mex('-v', '-outdir', progPath, ...
                    'CXXFLAGS=$CXXFLAGS -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include', ...
                    'LDFLAGS=$LDFLAGS -L/opt/homebrew/opt/libomp/lib', ...
                    'LINKLIBS=$LINKLIBS -lomp', ...
                    srcFile);
            end
        end
        %fprintf('done!\n');
    else
        fprintf('%s is up to date, skipping compilation.\n', ['semiCircle.' mexext]);
    end

end

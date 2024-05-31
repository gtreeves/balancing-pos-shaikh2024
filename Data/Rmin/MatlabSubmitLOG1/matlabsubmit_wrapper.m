% matlabsubmit wrapper function

%set the number of threads for the client
feature('NumThreads',1);
%set the the local profile and start the pool
lc=parcluster('local');
joblocation = strcat( pwd , '/MatlabSubmitLOG1');
lc.JobStorageLocation = joblocation;
if (~verLessThan('MATLAB','9.1'))
   lc.NumThreads=1;
end

lc.NumWorkers=45;
tamupool = parpool(lc,45);

if (verLessThan('MATLAB','9.1'))
spmd
   feature('NumThreads',1)
end
end
script_Smad_Rmin;
%close the local pool
delete(gcp('nocreate'));

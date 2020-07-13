function dump_spectra_init(dgngrp, currtime)
    FT = eltype(Settings.Q.data)
    bl = Settings.dg.balance_law
    mpicomm = Settings.mpicomm
    mpirank = MPI.Comm_rank(mpicomm)
    Q = Settings.Q
    spectrum, wavenumber = get_spectrum(Q, bl, dgngrp.interpol, dgngrp.nor, FT)
    if mpirank == 0
        dims = OrderedDict("k" => (wavenumber, Dict()))
        vars = OrderedDict("spectrum" => (("k",), FT, Dict()))

        dprefix = @sprintf(
            "%s_%s-%s",
            dgngrp.out_prefix,
            dgngrp.name,
            Settings.starttime,
        )
        dfilename = joinpath(Settings.output_dir, dprefix)
        init_data(dgngrp.writer, dfilename, dims, vars)
    end

    return nothing
end

function dump_spectra_collect(dgngrp, currtime)
    interpol = dgngrp.interpol
    mpicomm = Settings.mpicomm
    dg = Settings.dg
    Q = Settings.Q
    FT = eltype(Q.data)
    bl = dg.balance_law
    mpirank = MPI.Comm_rank(mpicomm)
    spectrum, wavenumber = get_spectrum(Q, bl, dgngrp.interpol, dgngrp.nor, FT)

    if mpirank == 0
        varvals = OrderedDict("spectrum" => spectrum)
        append_data(dgngrp.writer, varvals, currtime)
    end

    MPI.Barrier(mpicomm)
    return nothing
end

function dump_spectra_fini(dgngrp, currtime) end

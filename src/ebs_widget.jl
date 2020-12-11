using PyCall

widgets = PyNULL()
patches = PyNULL()
ticker = PyNULL()
gridspec = PyNULL()

function __init__()
    pygui(true)
    copy!(widgets, matplotlib.widgets)
    copy!(patches, matplotlib.patches)
    copy!(ticker, matplotlib.ticker)
    copy!(gridspec, matplotlib.gridspec)
end

function get_minmax(x::AbstractArray, r, xmin, xmax)
    _xmin, _xmax = extrema(x)
    if isnan(xmin) && !isnan(_xmin)
        xmin = _xmin
    else
        xmin = isnan(_xmin) ? xmin : min(xmin, _xmin)
    end
    if isnan(xmax) && !isnan(_xmax)
        xmax = _xmax
    else
        xmax = isnan(_xmax) ? xmax : max(xmax, _xmax)
    end
    return xmin - r, xmax + r
end

function get_minmax(x::AbstractArray, rx, y::AbstractArray, ry, xmin, xmax)
    xmin, xmax = get_minmax(x, rx, xmin, xmax)
    return get_minmax(y, ry, xmin, xmax)
end

function get_lims(xmin, xmax; pad=0.1)
    xpad = pad*(xmax - xmin)
    return xmin - xpad, xmax + xpad
end

get_lims(x::AbstractArray; kwargs...) = get_lims(extrema(x)...; kwargs...)

deg_formatter(x) = string("\$", Int(x), "^\\circ{}\$")
deg_formatter(x, pos) = deg_formatter(x)

function νlines(ax, ν, ω; lw=1)
    ν_deg = ustrip(u"°", ν)
    sep_pos = ax.axvline(ν_deg, color="tab:green", lw=lw)
    _sconj = ustrip(u"°", sconj(ω))
    supconj1 = ax.axvline(_sconj, ls="dashed", color="black", lw=lw)
    supconj2 = ax.axvline(_sconj - 360, ls="dashed", color="black", lw=lw)
    supconj3 = ax.axvline(_sconj + 360, ls="dashed", color="black", lw=lw)
    _iconj = ustrip(u"°", iconj(ω))
    infconj1 = ax.axvline(_iconj, ls="dotted", color="black", lw=lw)
    infconj2 = ax.axvline(_iconj - 360, ls="dotted", color="black", lw=lw)
    infconj3 = ax.axvline(_iconj + 360, ls="dotted", color="black", lw=lw)
    return sep_pos, (supconj1, supconj2, supconj3), (infconj1, infconj2, infconj3)
end

function νlines(sep_pos, supconj, infconj, ν, ω)
    ν_deg = ustrip(u"°", ν)
    sep_pos.set_xdata(ν_deg)
    _sconj = ustrip(u"°", sconj(ω))
    supconj[1].set_xdata(_sconj)
    supconj[2].set_xdata(_sconj - 360)
    supconj[3].set_xdata(_sconj + 360)
    _iconj = ustrip(u"°", iconj(ω))
    infconj[1].set_xdata(_iconj)
    infconj[2].set_xdata(_iconj - 360)
    infconj[3].set_xdata(_iconj + 360)
    return nothing
end

function νlines(xs::Vector{PyObject}, ν, ω)
    for x in xs
        νlines(x..., ν, ω)
    end
    return nothing
end

function get_xy_zneg(x, y, z)
    mask = @. sign(z) < 0
    return view(x, mask), view(y, mask)
end

function get_xy_zsky(x, y, z)
    mask = @. iszero(sign(z))
    return view(x, mask), view(y, mask)
end

function get_xy_zpos(x, y, z)
    mask = @. sign(z) > 0
    return view(x, mask), view(y, mask)
end

function ebs_widget(;
    figwidth = 10,
    figheight = 8,
    m1 = 2u"Msun",
    m2 = 1u"Msun",
    r1 = 2u"Rsun",
    r2 = 1u"Rsun",
    ν = 0u"rad",
    a = 0.1u"AU",
    e = 0,
    i = deg2rad(85)u"rad",
    ω = 0u"rad",
    Ω = 0u"rad",
    pcolor = "tab:orange",
    scolor = "tab:blue",
)

    r1_rsun = ustrip(u"Rsun", r1)
    r2_rsun = ustrip(u"Rsun", r2)
    a_rsun = ustrip(u"Rsun", a)
    ν_deg = ustrip(u"°", ν)
    ω_deg = ustrip(u"°", ω)

    fig = plt.figure(figsize=(figwidth, figheight), constrained_layout = true)
    gs0 = fig.add_gridspec(    # Define the outer gridspec
        nrows = 3,             # have sliders take up the bottom third so need three rows
        ncols = 2,             # plane of sky on left, separation plots on right
    )
    # separation plots
    gs1 = py"$gs0[0:2,1].subgridspec(2,1)"
    # sliders
    gs2 = py"$gs0[2,:].subgridspec(10,1)"

    ax_pos = fig.add_subplot(py"$gs0[0:2,0]")
    ax_pos.set_aspect("equal")
    ax_pos.set_anchor("W")
    ax_pos.spines["right"].set_position(("data", 0))
    ax_pos.spines["top"].set_position(("data", 0))

    ax_Δ² = fig.add_subplot(py"$gs1[0]")
    ax_dΔ² = fig.add_subplot(py"$gs1[1]", sharex = ax_Δ²)
    #ax_d²Δ² = fig.add_subplot(py"$gs1[2]", sharex = ax_Δ²)
    _xticks = 0:45:360
    ax_Δ².xaxis.set_ticks(_xticks)
    ax_Δ².xaxis.set_ticklabels(deg_formatter.(_xticks))

    axcolor = "lightgoldenrodyellow"
    ax_m1 = fig.add_subplot(py"$gs2[0]", fc = axcolor)
    ax_r1 = fig.add_subplot(py"$gs2[1]", fc = axcolor)
    ax_m2 = fig.add_subplot(py"$gs2[2]", fc = axcolor)
    ax_r2 = fig.add_subplot(py"$gs2[3]", fc = axcolor)
    ax_ν = fig.add_subplot(py"$gs2[4]", fc = axcolor)
    ax_a = fig.add_subplot(py"$gs2[5]", fc = axcolor)
    ax_e = fig.add_subplot(py"$gs2[6]", fc = axcolor)
    ax_i = fig.add_subplot(py"$gs2[7]", fc = axcolor)
    ax_ω = fig.add_subplot(py"$gs2[8]", fc = axcolor)
    ax_Ω = fig.add_subplot(py"$gs2[9]", fc = axcolor)
    sm1 = widgets.Slider(
        ax_m1, "\$M_1\$", 0.7, 10;
        valinit = ustrip(u"Msun", m1), closedmax=true, valfmt="%4.1f \$M_\\odot\$"
    )
    sm2 = widgets.Slider(
        ax_m2, "\$M_2\$", 0.7, 10;
        valinit = ustrip(u"Msun", m2), closedmax=true, valfmt="%4.1f \$M_\\odot\$"
    )
    sr1 = widgets.Slider(
        ax_r1, "\$R_1\$", 0, 10;
        valinit = ustrip(u"Rsun", r1), closedmax=false, valfmt="%4.1f \$R_\\odot\$"
    )
    sr2 = widgets.Slider(
        ax_r2, "\$R_2\$", 0, 10;
        valinit = ustrip(u"Rsun", r2), closedmax=false, valfmt="%4.1f \$R_\\odot\$"
    )
    sa = widgets.Slider(
        ax_a, "\$a\$", 0, 10;
        valinit = ustrip(u"AU", a), valstep=0.01, closedmax=false, valfmt="%5.2f AU"
    )
    se = widgets.Slider(ax_e, "\$e\$", 0, 1, valinit = e, valstep = 0.001, valfmt="%04.3f")
    sν = widgets.Slider(ax_ν, "\$\\nu\$", 0, 360; valinit = ustrip(u"°", ν), valfmt="%6.2f°")
    si = widgets.Slider(ax_i, "\$i\$", 0, 180; valinit = ustrip(u"°", i), valfmt="%6.2f°")
    sω = widgets.Slider(ax_ω, "\$\\omega\$", 0, 360; valinit = ustrip(u"°", ω), valfmt="%6.2f°")
    sΩ = widgets.Slider(ax_Ω, "\$\\Omega\$", 0, 360; valinit = ustrip(u"°", Ω), valfmt="%6.2f°")
    for s in [sr1, sr2, sν, sa, se, si, sω, sΩ]
        s.valtext.set_fontfamily("monospace")
    end

    nν = 1000
    νs = range(-ω, stop=(-ω + 2π*u"rad"), length=nν)
    
    xmin, xmax = Inf, -Inf
    ymin, ymax = Inf, -Inf

    kws = (
        lw = 1,
        marker="",
    )

    b = Binary(m1, r1, m2, r2, a; e=e, i=i, ω=ω, Ω=Ω)
    xyz1, xyz2 = pos(νs, b)

    x,y,z = ntuple(j -> ustrip.(u"Rsun", xyz1[j]), 3)
    vx, vy = get_xy_zneg(x, y, z)
    z1neg, = ax_pos.plot(vx, vy; ls = "dotted", zorder = -Inf, color = pcolor, kws...)
    vx, vy = get_xy_zsky(x, y, z)
    z1sky, = ax_pos.plot(vx, vy; ls = "dashed", zorder = 0, color = pcolor, kws...)
    vx, vy = get_xy_zpos(x, y, z)
    z1pos, = ax_pos.plot(vx, vy; ls="solid", zorder = eps(), color = pcolor, kws...)
    xmin, xmax = get_minmax(x, r1_rsun, xmin, xmax)
    ymin, ymax = get_minmax(y, r1_rsun, ymin, ymax)

    x,y,z = ntuple(j -> ustrip.(u"Rsun", xyz2[j]), 3)
    vx, vy = get_xy_zneg(x, y, z)
    z2neg, = ax_pos.plot(vx, vy; ls = "dotted", zorder = -Inf, color = scolor, kws...)
    vx, vy = get_xy_zsky(x, y, z)
    z2sky, = ax_pos.plot(vx, vy; ls = "dashed", zorder = 0, color = scolor, kws...)
    vx, vy = get_xy_zpos(x, y, z)
    z2pos, = ax_pos.plot(vx, vy; ls = "solid", zorder = eps(), color = scolor, kws...)
    xmin, xmax = get_minmax(x, r2_rsun, xmin, xmax)
    ymin, ymax = get_minmax(y, r2_rsun, ymin, ymax)

    xmin, xmax = get_lims(xmin, xmax)
    ax_pos.set_xlim(xmin, xmax)
    ax_pos.set_ylim(get_lims(ymin, ymax)...)
    ax_pos.add_patch(
        patches.Rectangle((0,0), 1, 1;
            transform = ax_pos.transAxes,
            zorder = 0,
            color = "black",
            alpha = 0.2
        )
    )

    xyz1, xyz2 = pos(ν, b)
    x,y,z = ustrip.(u"Rsun", xyz1)
    pri = ax_pos.add_patch(
        patches.Circle(
            (x,y),
            r1_rsun,
            facecolor = pcolor,
            edgecolor = "none",
            zorder = z
        )
    )
    x,y,z = ustrip.(u"Rsun", xyz2)
    sec = ax_pos.add_patch(
        patches.Circle(
            (x,y),
            r2_rsun,
            facecolor = scolor,
            edgecolor = "none",
            zorder = z
        )
    )

    νmin, νmax = -π*u"rad"/6, 13π*u"rad"/6
    νs_sep = range(νmin, stop=νmax, length=1000)
    ds_sep = ustrip.(u"°", νs_sep)
    νmin_deg, νmax_deg = ustrip(u"°", νmin), ustrip(u"°", νmax)
    
    ax_Δ².set_ylabel("\$\\Delta^2 \\ \\left[\\mathrm{R}_\\odot^2\\right]\$")
    fy_rsun = ustrip.(u"Rsun^2", Δ²(νs_sep, b))
    fplot, = ax_Δ².plot(ds_sep, fy_rsun; color="black", kws...)
    fr1r2 = ax_Δ².axhline((r1_rsun^2 + r2_rsun^2)/a_rsun^2, color="black", ls="dotted")
    _, fy_max = get_lims(0, maximum(fy_rsun))
    ax_Δ².set_xlim(νmin_deg, νmax_deg)
    ax_Δ².set_ylim(0, fy_max)

    ax_dΔ².set_ylabel(
        "\$"*
        "\\frac{\\mathrm{d}}{\\mathrm{d}\\nu}"*
        "\\left(\\Delta^2\\right)"*
        "\$"
    )
    gy_rsun = ustrip.(u"Rsun^2", dΔ²_dν(νs_sep, b))
    gplot, = ax_dΔ².plot(ds_sep, gy_rsun; color="black", kws...)
    ax_dΔ².set_ylim(get_lims(gy_rsun)...)
    ax_dΔ².spines["bottom"].set_position(("data", 0))

    #ax_d²Δ².set_xlabel("\$\\nu\$")
    #ax_d²Δ².set_ylabel(
    #    "\$"*
    #    "\\frac{\\mathrm{d}^2}{\\mathrm{d}^2\\nu}"*
    #    "\\left(\\Delta^2\\right)"*
    #    "\$"
    #)
    #hy_rsun = ustrip.(u"Rsun^2", d²Δ²_dν²(νs_sep, b))
    #hplot, = ax_d²Δ².plot(ds_sep, hy_rsun; color="black", kws...)
    ##ax_d²Δ².set_ylim(get_lims(hy)...)
    #ax_d²Δ².spines["bottom"].set_position(("data", 0))
    
    # current position (in terms of true anomaly, supconj and infconj)
    _νlines = Vector{PyObject}(undef, 2)
    _νspans = Vector{PyObject}(undef, 4)
    for (j,ax) in enumerate([ax_Δ², ax_dΔ²])
        _νlines[j] = νlines(ax, ν, ω)
        _νspans[j] = ax.axvspan(νmin_deg, 0; color="black", alpha=0.1)
        _νspans[j+2] = ax.axvspan(360, νmax_deg; color="black", alpha=0.1)
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
    end

    ν₁, ν₂ = get_eclipses_ν(b)
    ν₁_deg, ν₂_deg = ustrip.(u"°", (ν₁, ν₂))
    kws = (linewidth = 2, alpha = 0.5, zorder = -1)
    Δ²pecl = ax_Δ².axvline(ν₁_deg; color=pcolor, kws...)
    Δ²secl = ax_Δ².axvline(ν₂_deg; color=scolor, kws...)
    dΔ²pecl = ax_dΔ².axvline(ν₁_deg; color=pcolor, kws...)
    dΔ²secl = ax_dΔ².axvline(ν₂_deg; color=scolor, kws...)
    #d²Δ²pecl = ax_d²Δ².axvline(νdeg₁; color=pcolor, kws...)
    #d²Δ²secl = ax_d²Δ².axvline(νdeg₂; color=scolor, kws...)

    function update_νonly(val)
        m1_msun = sm1.val
        m1 = m1_msun*u"Msun"
        m2_msun = sm2.val
        m2 = m2_msun*u"Msun"
        r1_rsun = sr1.val
        r1 = r1_rsun*u"Rsun"
        r2_rsun = sr2.val
        r2 = r2_rsun*u"Rsun"
        a_au = sa.val
        a = a_au*u"AU"

        e = se.val

        ν_deg = sν.val
        ν = ν_deg*u"°"
        i_deg = si.val
        i = i_deg*u"°"
        ω_deg = sω.val
        ω = ω_deg*u"°"
        Ω_deg = sΩ.val
        Ω = Ω_deg*u"°"

        b = Binary(m1, r1, m2, r2, a; e=e, i=i, ω=ω, Ω=Ω)
        
        xyz1, xyz2 = pos(ν, b)
        x,y,z = ustrip.(u"Rsun", xyz1)
        pri.set_center((x, y))
        pri.set_zorder(z)
        pri.set_radius(r1_rsun)
        x,y,z = ustrip.(u"Rsun", xyz2)
        sec.set_center((x, y))
        sec.set_zorder(z)
        sec.set_radius(r2_rsun)

        νlines(_νlines, ν, ω)
    end

    function update(val)
        m1_msun = sm1.val
        m1 = m1_msun*u"Msun"
        m2_msun = sm2.val
        m2 = m2_msun*u"Msun"
        r1_rsun = sr1.val
        r1 = r1_rsun*u"Rsun"
        r2_rsun = sr2.val
        r2 = r2_rsun*u"Rsun"
        a_au = sa.val
        a = a_au*u"AU"

        e = se.val

        ν_deg = sν.val
        ν = ν_deg*u"°"
        i_deg = si.val
        i = i_deg*u"°"
        ω_deg = sω.val
        ω = ω_deg*u"°"
        Ω_deg = sΩ.val
        Ω = Ω_deg*u"°"

        b = Binary(m1, r1, m2, r2, a; e=e, i=i, ω=ω, Ω=Ω)
       
        xmin, xmax = Inf, -Inf
        ymin, ymax = Inf, -Inf

        νs = range(-ω, stop=(-ω + 360u"°"), length=nν)

        xyz1, xyz2 = pos(νs, b)
        x,y,z = ntuple(j -> ustrip.(u"Rsun", xyz1[j]), 3)
        vx, vy = get_xy_zneg(x, y, z)
        z1neg.set_data(vx, vy)
        vx, vy = get_xy_zsky(x, y, z)
        z1sky.set_data(vx, vy)
        vx, vy = get_xy_zpos(x, y, z)
        z1pos.set_data(vx, vy)
        xmin, xmax = get_minmax(x, r1_rsun, xmin, xmax)
        ymin, ymax = get_minmax(y, r1_rsun, ymin, ymax)
        x,y,z = ntuple(j -> ustrip.(u"Rsun", xyz2[j]), 3)
        vx, vy = get_xy_zneg(x, y, z)
        z2neg.set_data(vx, vy)
        vx, vy = get_xy_zsky(x, y, z)
        z2sky.set_data(vx, vy)
        vx, vy = get_xy_zpos(x, y, z)
        z2pos.set_data(vx, vy)
        xmin, xmax = get_minmax(x, r2_rsun, xmin, xmax)
        ymin, ymax = get_minmax(y, r2_rsun, ymin, ymax)

        update_νonly(nothing)

        ax_pos.set_xlim(get_lims(xmin, xmax)...)
        ax_pos.set_ylim(get_lims(ymin, ymax)...)


        fy_rsun = ustrip.(u"Rsun^2", Δ²(νs_sep, b))
        fplot.set_ydata(fy_rsun)
        _, fy_max = get_lims(0, maximum(fy_rsun))
        ax_Δ².set_ylim(0, fy_max)

        gy_rsun = ustrip.(u"Rsun^2", dΔ²_dν(νs_sep, b))
        gplot.set_ydata(gy_rsun)
        ax_dΔ².set_ylim(get_lims(gy_rsun)...)

        ν₁, ν₂ = get_eclipses_ν(b)
        ν₁_deg, ν₂_deg = ustrip.(u"°", (ν₁, ν₂))
        Δ²pecl.set_xdata(ν₁_deg)
        dΔ²pecl.set_xdata(ν₁_deg)
        #d²Δ²pecl.set_xdata(νdeg)
        Δ²secl.set_xdata(ν₂_deg)
        dΔ²secl.set_xdata(ν₂_deg)
        #d²Δ²secl.set_xdata(νdeg)
        νlines(_νlines, ν, ω)
    end
    #
    sm1.on_changed(update)
    sm2.on_changed(update)
    sr1.on_changed(update)
    sr2.on_changed(update)
    sν.on_changed(update_νonly)
    sa.on_changed(update)
    se.on_changed(update)
    si.on_changed(update)
    sω.on_changed(update)
    sΩ.on_changed(update)
    
    nothing
end

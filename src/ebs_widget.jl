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

function get_lims(xmin, xmax; pad=0.1)
    xpad = pad*(xmax - xmin)
    return xmin - xpad, xmax + xpad
end

get_lims(x::AbstractArray; kwargs...) = get_lims(extrema(x)...; kwargs...)

#function get_lims(x::AbstractArray, r; kwargs...)
#    xmin, xmax = get_lims(x)
#    return xmin - r, xmax + r
#end

deg_formatter(x) = string("\$", Int(x), "^\\circ{}\$")
deg_formatter(x, pos) = deg_formatter(x)

function νlines(ax, ν, ω; lw=1)
    sep_pos = ax.axvline(rad2deg(ν), color="tab:green", lw=lw)
    _sconj = rad2deg(sconj(ω))
    supconj1 = ax.axvline(_sconj, ls="dashed", color="black", lw=lw)
    supconj2 = ax.axvline(_sconj - 360, ls="dashed", color="black", lw=lw)
    supconj3 = ax.axvline(_sconj + 360, ls="dashed", color="black", lw=lw)
    _iconj = rad2deg(iconj(ω))
    infconj1 = ax.axvline(_iconj, ls="dotted", color="black", lw=lw)
    infconj2 = ax.axvline(_iconj - 360, ls="dotted", color="black", lw=lw)
    infconj3 = ax.axvline(_iconj + 360, ls="dotted", color="black", lw=lw)
    return sep_pos, (supconj1, supconj2, supconj3), (infconj1, infconj2, infconj3)
end

function νlines(sep_pos, supconj, infconj, ν, ω)
    sep_pos.set_xdata(rad2deg(ν))
    _sconj = rad2deg(sconj(ω))
    supconj[1].set_xdata(_sconj)
    supconj[2].set_xdata(_sconj - 360)
    supconj[3].set_xdata(_sconj + 360)
    _iconj = rad2deg(iconj(ω))
    infconj[1].set_xdata(_iconj)
    infconj[2].set_xdata(_iconj - 360)
    infconj[3].set_xdata(_iconj + 360)
    nothing
end

function νlines(xs::Vector{PyObject}, ν, ω)
    for x in xs
        νlines(x..., ν, ω)
    end
    nothing
end

#function ebs_test()
#    m1 = 1*u"Msun"
#    m2 = 1*u"Msun"
#    r1 = 1*u"Rsun"
#    r2 = 1*u"Rsun"
#    ν = 0*u"°"
#    a = 0.1*u"AU"
#    e = 0
#    i = 85*u"°"
#    ω = 0*u"°"
#    Ω = 0*u"°"
#
#    s1 = Star(m1, r1)
#    s2 = Star(m2, r2)
#    orb = Orbit(s1, s2, a; e=e, i=i, ω=ω, Ω=Ω)
#    bin = Binary(s1, s2, orb)
#    EclipsingBinary(bin)
#end

reduced_mass(M1, M2) = M1*M2/(M1 + M2)

function pos1(ν, m1, m2, a, e, i, ω, Ω)
    a = a*reduced_mass(m1, m2)/m1
    return pos(ν, a, e, i, ω, Ω)
end

function pos2(ν, m1, m2, a, e, i, ω, Ω)
    a = -a*reduced_mass(m1, m2)/m2
    return pos(ν, a, e, i, ω, Ω)
end

function get_xy_zneg(x, y, z)
    mask = @. sign(z) < 0
    return view(x, mask), view(y, mask)
end

function get_xy_zpos(x, y, z)
    mask = @. sign(z) > 0
    return view(x, mask), view(y, mask)
end

function get_xy_zsky(x, y, z)
    mask = @. iszero(sign(z))
    return view(x, mask), view(y, mask)
end

function ebs_widget(;
    figwidth = 10,
    figheight = 8,
    m1 = 2,
    m2 = 1,
    r1 = 2,
    r2 = 1,
    ν = 0,
    a = 0.1,
    e = 0,
    i = deg2rad(85),
    ω = 0,
    Ω = 0,
    pcolor = "tab:orange",
    scolor = "tab:blue",
)
    fig = plt.figure(figsize=(figwidth, figheight), constrained_layout = true)
    gs0 = fig.add_gridspec(    # Define the outer gridspec
        nrows = 3,             # have sliders take up the bottom third so need three rows
        ncols = 2,             # plane of sky on left, separation plots on right
    )
    # separation plots
    gs1 = py"$gs0[0:2,1].subgridspec(3,1)"
    # sliders
    gs2 = py"$gs0[2,:].subgridspec(10,1)"

    ax_pos = fig.add_subplot(py"$gs0[0:2,0]")
    ax_pos.set_aspect("equal")
    #ax_pos.set_anchor("W")
    ax_pos.spines["right"].set_position(("data", 0))
    ax_pos.spines["top"].set_position(("data", 0))

    ax_Δ² = fig.add_subplot(py"$gs1[0]")
    ax_dΔ² = fig.add_subplot(py"$gs1[1]", sharex = ax_Δ²)
    ax_d²Δ² = fig.add_subplot(py"$gs1[2]", sharex = ax_Δ²)
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
    sm1 = widgets.Slider(ax_m1, "\$M_1\$", 0.7, 10;
        valinit = m1, closedmax=true, valfmt="%4.1f \$M_\\odot\$"
    )
    sr1 = widgets.Slider(ax_r1, "\$R_1\$", 0, 10;
        valinit = r1, closedmax=false, valfmt="%4.1f \$R_\\odot\$"
    )
    sm2 = widgets.Slider(ax_m2, "\$M_2\$", 0.7, 10;
        valinit = m2, closedmax=true, valfmt="%4.1f \$M_\\odot\$"
    )
    sr2 = widgets.Slider(ax_r2, "\$R_2\$", 0, 10;
        valinit = r2, closedmax=false, valfmt="%4.1f \$R_\\odot\$"
    )
    sν = widgets.Slider(ax_ν, "\$\\nu\$", 0, 360, valinit = rad2deg(ν), valfmt="%6.2f°")
    sa = widgets.Slider(ax_a, "\$a\$", 0, 100;
        valinit = a, valstep=0.1, closedmax=false, valfmt="%4.1f AU"
    )
    a = ustrip(u"Rsun", a*u"AU")

    se = widgets.Slider(ax_e, "\$e\$", 0, 1, valinit = e, valstep = 0.001, valfmt="%04.3f")
    si = widgets.Slider(ax_i, "\$i\$", 0, 180, valinit = rad2deg(i), valfmt="%6.2f°")
    sω = widgets.Slider(ax_ω, "\$\\omega\$", 0, 360, valinit = rad2deg(ω), valfmt="%6.2f°")
    sΩ = widgets.Slider(ax_Ω, "\$\\Omega\$", 0, 360, valinit = rad2deg(Ω), valfmt="%6.2f°")
    for s in [sr1, sr2, sν, sa, se, si, sω, sΩ]
        s.valtext.set_fontfamily("monospace")
    end

    nν = 1000
    νs = range(-ω, stop=(-ω + 2π), length=nν)
    
    xmin, xmax = Inf, -Inf
    ymin, ymax = Inf, -Inf

    kws = (lw = 1, marker="")
    x, y, z = pos1(νs, m1, m2, a, e, i, ω, Ω)
    vx, vy = get_xy_zneg(x, y, z)
    z1neg, = ax_pos.plot(vx, vy; ls = "dotted", zorder = -10, color = pcolor, kws...)
    vx, vy = get_xy_zsky(x, y, z)
    z1sky, = ax_pos.plot(vx, vy; ls = "dashed", zorder = 0, color = pcolor, kws...)
    vx, vy = get_xy_zpos(x, y, z)
    z1pos, = ax_pos.plot(vx, vy; ls="solid", zorder=eps(), color = pcolor, kws...)
    xmin, xmax = get_minmax(x, r1, xmin, xmax)
    ymin, ymax = get_minmax(y, r1, ymin, ymax)

    x, y, z = pos2(νs, m1, m2, a, e, i, ω, Ω)
    vx, vy = get_xy_zneg(x, y, z)
    z2neg, = ax_pos.plot(vx, vy; ls="dotted", zorder=-10, color = scolor, kws...)
    vx, vy = get_xy_zsky(x, y, z)
    z2sky, = ax_pos.plot(vx, vy; ls="dashed", zorder=0, color = scolor, kws...)
    vx, vy = get_xy_zpos(x, y, z)
    z2pos, = ax_pos.plot(vx, vy; ls="solid", zorder=eps(), color = scolor, kws...)
    xmin, xmax = get_minmax(x, r2, xmin, xmax)
    ymin, ymax = get_minmax(y, r2, ymin, ymax)

    ax_pos.set_xlim(get_lims(xmin, xmax)...)
    ax_pos.set_ylim(get_lims(ymin, ymax)...)

    x, y, z = pos1(ν, m1, m2, a, e, i, ω, Ω)
    pri = ax_pos.add_patch(
        patches.Circle(
            (x,y),
            r1,
            facecolor = pcolor,
            edgecolor = "none",
            zorder = z
        )
    )
    x, y, z = pos2(ν, m1, m2, a, e, i, ω, Ω)
    sec = ax_pos.add_patch(
        patches.Circle(
            (x,y),
            r2,
            facecolor = scolor,
            edgecolor = "none",
            zorder = z
        )
    )


    νmin, νmax = -π/6, 13π/6
    νs_sep = range(νmin, stop=νmax, length=1000)
    ds_sep = rad2deg.(νs_sep)
    νdegmin, νdegmax = rad2deg(νmin), rad2deg(νmax)
    
    fy = Δ²(νs_sep, a, e, i, ω)
    fplot, = ax_Δ².plot(ds_sep, fy; color="black", kws...)
    fr1r2 = ax_Δ².axhline((r1^2 + r2^2)/a^2, color="black", ls="dotted")
    #ax_Δ².set_xlim(extrema(ds_sep)...)
    #ax_Δ².set_ylim(0, 1.1*maximum(fy))
    ax_Δ².set_ylabel("\$\\Delta^2 \\ [\\mathrm{a}^2]\$")

    ax_dΔ².set_ylabel(
        "\$"*
        "\\frac{\\mathrm{d}}{\\mathrm{d}\\nu}"*
        "\\left(\\Delta^2\\right)"*
        "\$"
    )
    gy = dΔ²_dν(νs_sep, a, e, i, ω)
    gplot, = ax_dΔ².plot(ds_sep, gy; color="black", kws...)
    ax_dΔ².set_ylim(get_lims(gy)...)
    ax_dΔ².spines["bottom"].set_position(("data", 0))

    ax_d²Δ².set_xlabel("\$\\nu\$")
    ax_d²Δ².set_ylabel(
        "\$"*
        "\\frac{\\mathrm{d}^2}{\\mathrm{d}^2\\nu}"*
        "\\left(\\Delta^2\\right)"*
        "\$"
    )
    hy = d²Δ²_dν²(νs_sep, a, e, i, ω)
    hplot, = ax_d²Δ².plot(ds_sep, hy; color="black", kws...)
    ax_d²Δ².set_ylim(get_lims(hy)...)
    ax_d²Δ².spines["bottom"].set_position(("data", 0))
    
    # current position (in terms of true anomaly, supconj and infconj)
    _νlines = Vector{PyObject}(undef, 3)
    _νspans = Vector{PyObject}(undef, 6)
    for (j,ax) in enumerate([ax_Δ², ax_dΔ², ax_d²Δ²])
        _νlines[j] = νlines(ax, ν, ω)
        _νspans[j] = ax.axvspan(νdegmin, 0; color="black", alpha=0.1)
        _νspans[j+3] = ax.axvspan(360, νdegmax; color="black", alpha=0.1)
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.set_xlim(νdegmin, νdegmax)
    end

    ν₁, ν₂ = get_eclipses_ν(r1, r2, a, e, i, ω)
    νdeg₁, νdeg₂ = ustrip.(u"°", (ν₁, ν₂))
    kws = (linewidth = 2, alpha = 0.5, zorder = -1)
    Δ²pecl = ax_Δ².axvline(νdeg₁; color=pcolor, kws...)
    Δ²secl = ax_Δ².axvline(νdeg₂; color=scolor, kws...)
    dΔ²pecl = ax_dΔ².axvline(νdeg₁; color=pcolor, kws...)
    dΔ²secl = ax_dΔ².axvline(νdeg₂; color=scolor, kws...)
    d²Δ²pecl = ax_d²Δ².axvline(νdeg₁; color=pcolor, kws...)
    d²Δ²secl = ax_d²Δ².axvline(νdeg₂; color=scolor, kws...)


    #x, y, z = pos(ν₁, a, e, i, ω, Ω)
    #_x, _y, _z = pos1((x, y, z), m1, m2)
    #ls = _z > 0 ? "solid" : "dotted"
    #pecl1 = ax_pos.add_patch(
    #    patches.Circle(
    #        (_x, _y),
    #        r1;
    #        facecolor = "none",
    #        edgecolor = pcolor,
    #        linestyle = ls,
    #        zorder = _z
    #    )
    #)
    #_x, _y, _z = pos2((x, y, z), m1, m2)
    #ls = _z > 0 ? "solid" : "dotted"
    #secl1 = ax_pos.add_patch(
    #    patches.Circle(
    #        (_x, _y),
    #        r2;
    #        facecolor = "none",
    #        edgecolor = scolor,
    #        linestyle = ls,
    #        zorder = _z
    #    )
    #)

    #x, y, z = pos(ν₂, a, e, i, ω, Ω)
    #_x, _y, _z = pos1((x, y, z), m1, m2)
    #ls = _z > 0 ? "solid" : "dotted"
    #pecl2 = ax_pos.add_patch(
    #    patches.Circle(
    #        (_x, _y),
    #        r1;
    #        facecolor = "none",
    #        edgecolor = pcolor,
    #        linestyle = ls,
    #        zorder = _z
    #    )
    #)
    #_x, _y, _z = pos2((x, y, z), m1, m2)
    #ls = _z > 0 ? "solid" : "dotted"
    #secl2 = ax_pos.add_patch(
    #    patches.Circle(
    #        (_x, _y),
    #        r2;
    #        facecolor = "none",
    #        edgecolor = scolor,
    #        linestyle = ls,
    #        zorder = _z
    #    )
    #)

    function update_νonly(val)
        m1 = sm1.val
        m2 = sm2.val
        r1 = sr1.val
        r2 = sr2.val
        ν = deg2rad(sν.val)
        a = ustrip(u"Rsun", sa.val*u"AU")
        e = se.val
        i = deg2rad(si.val)
        ω = deg2rad(sω.val)
        Ω = deg2rad(sΩ.val)
        
        x, y, z = pos1(ν, m1, m2, a, e, i, ω, Ω)
        pri.set_center((x, y))
        pri.set_zorder(z)
        pri.set_radius(r1)
        x, y, z = pos2(ν, m1, m2, a, e, i, ω, Ω)
        sec.set_center((x, y))
        sec.set_zorder(z)
        sec.set_radius(r2)

        νlines(_νlines, ν, ω)
    end

    function update(val)
        m1 = sm1.val
        m2 = sm2.val
        r1 = sr1.val
        r2 = sr2.val
        ν = deg2rad(sν.val)
        a = ustrip(u"Rsun", sa.val*u"AU")
        e = se.val
        i = deg2rad(si.val)
        ω = deg2rad(sω.val)
        Ω = deg2rad(sΩ.val)
        
        xmin, xmax = Inf, -Inf
        ymin, ymax = Inf, -Inf

        νs = range(-ω, stop=(-ω + 2π), length=nν)

        x, y, z = pos1(νs, m1, m2, a, e, i, ω, Ω)
        xmin, xmax = get_minmax(x, r1, xmin, xmax)
        ymin, ymax = get_minmax(y, r1, ymin, ymax)
        vx, vy = get_xy_zneg(x, y, z)
        z1neg.set_data(vx, vy)
        vx, vy = get_xy_zsky(x, y, z)
        z1sky.set_data(vx, vy)
        vx, vy = get_xy_zpos(x, y, z)
        z1pos.set_data(vx, vy)

        x, y, z = pos2(νs, m1, m2, a, e, i, ω, Ω)
        xmin, xmax = get_minmax(x, r2, xmin, xmax)
        ymin, ymax = get_minmax(y, r2, ymin, ymax)
        vx, vy = get_xy_zneg(x, y, z)
        z2neg.set_data(vx, vy)
        vx, vy = get_xy_zsky(x, y, z)
        z2sky.set_data(vx, vy)
        vx, vy = get_xy_zpos(x, y, z)
        z2pos.set_data(vx, vy)

        update_νonly(nothing)

        ax_pos.set_xlim(get_lims(xmin, xmax)...)
        ax_pos.set_ylim(get_lims(ymin, ymax)...)

        ν₁, ν₂ = get_eclipses_ν(r1, r2, a, e, i, ω)
        νdeg = ustrip(u"°", ν₁)
        Δ²pecl.set_xdata(νdeg)
        dΔ²pecl.set_xdata(νdeg)
        d²Δ²pecl.set_xdata(νdeg)
        νdeg = ustrip(u"°", ν₂)
        Δ²secl.set_xdata(νdeg)
        dΔ²secl.set_xdata(νdeg)
        d²Δ²secl.set_xdata(νdeg)
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

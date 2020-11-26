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

function get_minmax(x::AbstractArray, xmin, xmax)
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
    return xmin, xmax
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

function ebs_test()
    m1 = 1*u"Msun"
    m2 = 1*u"Msun"
    r1 = 1*u"Rsun"
    r2 = 1*u"Rsun"
    ν = 0*u"°"
    a = 0.1*u"AU"
    e = 0
    i = 85*u"°"
    ω = 0*u"°"
    Ω = 0*u"°"

    s1 = Star(m1, r1)
    s2 = Star(m2, r2)
    orb = Orbit(s1, s2, a; e=e, i=i, ω=ω, Ω=Ω)
    bin = Binary(s1, s2, orb)
    EclipsingBinary(bin)
end

function ebs_widget(;
    figwidth = 10,
    figheight = 8,
    r1 = 1,
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
    gs2 = py"$gs0[2,:].subgridspec(8,1)"

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
    ax_r1 = fig.add_subplot(py"$gs2[0]", fc = axcolor)
    ax_r2 = fig.add_subplot(py"$gs2[1]", fc = axcolor)
    ax_ν = fig.add_subplot(py"$gs2[2]", fc = axcolor)
    ax_a = fig.add_subplot(py"$gs2[3]", fc = axcolor)
    ax_e = fig.add_subplot(py"$gs2[4]", fc = axcolor)
    ax_i = fig.add_subplot(py"$gs2[5]", fc = axcolor)
    ax_ω = fig.add_subplot(py"$gs2[6]", fc = axcolor)
    ax_Ω = fig.add_subplot(py"$gs2[7]", fc = axcolor)
    sr1 = widgets.Slider(ax_r1, "\$R_1 \\ [a]\$", 0, 10;
        valinit = r1, closedmax=false, valfmt="%4.1f \$R_\\odot\$"
    )
    sr2 = widgets.Slider(ax_r2, "\$R_2 \\ [a]\$", 0, 10;
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
    x, y, z = pos(νs, a, e, i, ω, Ω)
    #_x = ustrip.(u"Rsun", x)
    #_y = ustrip.(u"Rsun", y)
    
    #ax1.set_xlabel("\$x \\ [\\mathrm{a}]\$", ha="right", va="bottom")
    #ax1.set_ylabel("\$y \\ [\\mathrm{a}]\$", ha="right", va="top")

    kws = (color = "black", lw = 1, marker="")
    mask = @. sign(z) < 0
    zneg, = ax_pos.plot(view(x, mask), view(y, mask); ls="dotted", zorder=-10, kws...)
    mask = @. iszero(sign(z))
    zsky, = ax_pos.plot(view(x, mask), view(y, mask); ls="dashed", zorder=0, kws...)
    mask = @. sign(z) > 0
    zpos, = ax_pos.plot(view(x, mask), view(y, mask); ls="solid", zorder=eps(), kws...)

    x, y, z = pos(ν, a, e, i, ω, Ω)
    
    pri = ax_pos.add_patch(
        patches.Circle(
            (0,0),
            r1,
            facecolor = pcolor,
            edgecolor = "none",
            zorder = 0
        )
    )
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
    fplot, = ax_Δ².plot(ds_sep, fy; kws...)
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
    gplot, = ax_dΔ².plot(ds_sep, gy; kws...)
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
    hplot, = ax_d²Δ².plot(ds_sep, hy; kws...)
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

    ν₁, ν₂ = get_eclipses(r1, r2, a, e, i, ω)
    x, y, z = pos(ν₁, a, e, i, ω, Ω)
    pecl1 = ax_pos.add_patch(
        patches.Circle(
            (x, y),
            r2;
            facecolor = scolor,
            edgecolor = "none",
            alpha = 0.5,
            zorder = z
        )
    )
    x, y, z = pos(ν₂, a, e, i, ω, Ω)
    secl1 = ax_pos.add_patch(
        patches.Circle(
            (x,y),
            r2;
            facecolor = scolor,
            edgecolor = "none",
            alpha = 0.5,
            zorder = z
        )
    )

    ν₁, ν₂ = rad2deg.((ν₁, ν₂))
    kws = (linewidth = 2, alpha = 0.5, zorder = -1)
    pecl2 = ax_Δ².axvline(ν₁; color=pcolor, kws...)
    pecl3 = ax_dΔ².axvline(ν₁; color=pcolor, kws...)
    pecl4 = ax_d²Δ².axvline(ν₁; color=pcolor, kws...)
    secl2 = ax_Δ².axvline(ν₂; color=scolor, kws...)
    secl3 = ax_dΔ².axvline(ν₂; color=scolor, kws...)
    secl4 = ax_d²Δ².axvline(ν₂; color=scolor, kws...)

    function update(val)
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
        if iszero(i)
            νs = range(0, stop=2π, length=nν)
            x, y, _ = pos(νs, a, e, i, ω, Ω)
            xmin, xmax = get_minmax(x, xmin, xmax)
            ymin, ymax = get_minmax(y, ymin, ymax)
            zsky.set_data(x, y)
            zneg.set_data([], [])
            zpos.set_data([], [])
            pecl1.set_center((NaN,NaN))
            secl1.set_center((NaN,NaN))
            pecl2.set_xdata(NaN)
            secl2.set_xdata(NaN)
            pecl3.set_xdata(NaN)
            secl3.set_xdata(NaN)
            pecl4.set_xdata(NaN)
            secl4.set_xdata(NaN)
            #ax_pos.set_xlim(get_lims(x, r2)...)
            #ax_pos.set_ylim(get_lims(y, r2)...)
        else
            ν₁, ν₂ = get_eclipses(r1, r2, a, e, i, ω)

            x, y, z = pos(ν₁, a, e, i, ω, Ω)
            pecl1.set_center((x,y))
            pecl1.set_zorder(z)
            pecl1.set_radius(r2)
            xmin = min(xmin, x - r2)
            xmax = max(xmax, x + r2)
            ymin = min(ymin, y - r2)
            ymax = max(ymax, y + r2)
            νd₁ = rad2deg(ν₁)
            pecl2.set_xdata(νd₁)
            pecl3.set_xdata(νd₁)
            pecl4.set_xdata(νd₁)
            
            x, y, z = pos(ν₂, a, e, i, ω, Ω)
            secl1.set_center((x,y))
            secl1.set_zorder(z)
            secl1.set_radius(r2)
            xmin = min(xmin, x - r2)
            xmax = max(xmax, x + r2)
            ymin = min(ymin, y - r2)
            ymax = max(ymax, y + r2)
            νd₂ = rad2deg(ν₂)
            secl2.set_xdata(νd₂)
            secl3.set_xdata(νd₂)
            secl4.set_xdata(νd₂)
            
            if i < π
                zsky.set_data([], [])
                νs = range(-ω, stop=(π-ω), length=nν)
                x, y, z = pos(νs, a, e, i, ω, Ω)
                zneg.set_data(x, y)
                xmin, xmax = get_minmax(x, xmin, xmax)
                ymin, ymax = get_minmax(y, ymin, ymax)

                νs = range(π-ω, stop=(2π-ω), length=nν)
                x, y, _ = pos(νs, a, e, i, ω, Ω)
                zpos.set_data(x, y)
                xmin, xmax = get_minmax(x, xmin, xmax)
                ymin, ymax = get_minmax(y, ymin, ymax)
            else
                zsky.set_data([], [])
                νs = range(-ω, stop=(π-ω), length=nν)
                x, y, z = pos(νs, a, e, i, ω, Ω)
                zpos.set_data(x, y)
                xmin, xmax = get_minmax(x, xmin, xmax)
                ymin, ymax = get_minmax(y, ymin, ymax)

                νs = range(π-ω, stop=(2π-ω), length=nν)
                x, y, _ = pos(νs, a, e, i, ω, Ω)
                zneg.set_data(x, y)
                xmin, xmax = get_minmax(x, xmin, xmax)
                ymin, ymax = get_minmax(y, ymin, ymax)
            end
        end
        x, y, z = pos(ν, a, e, i, ω, Ω)
        sec.set_center((x, y))
        sec.set_zorder(z)
        pri.set_radius(r1)
        sec.set_radius(r2)

        xmin = min(xmin, -r1, x - r2)
        xmax = max(xmax, r1, x + r2)
        ymin = min(ymin, -r1, y - r2)
        ymax = max(ymax, r1, y + r2)

        ax_pos.set_xlim(get_lims(xmin, xmax))
        ax_pos.set_ylim(get_lims(ymin, ymax))
        
        fy = Δ²(νs_sep, a, e, i, ω)
        fplot.set_data(ds_sep, fy)
        ax_Δ².set_ylim(get_lims(fy)...)
        
        gy = dΔ²_dν(νs_sep, a, e, i, ω)
        gplot.set_data(ds_sep, gy)
        ax_dΔ².set_ylim(get_lims(gy)...)
        
        hy = d²Δ²_dν²(νs_sep, a, e, i, ω)
        hplot.set_data(ds_sep, hy)
        ax_d²Δ².set_ylim(get_lims(hy)...)

        νlines(_νlines, ν, ω)
        #if iszero(i) && iszero(e)
        #    ν_pe, ν_se = NaN, NaN
        #else
        #    ν_pe, ν_se = get_eclipses(e, i, ω)
        #end
        #ν_pe, ν_se = get_minimums(e, i, ω)
        #eνlines(eclips1..., ν_pe, ν_se)
        #eνlines(eclips2..., ν_pe, ν_se)
        #eνlines(eclips3..., ν_pe, ν_se)
        fig.canvas.draw_idle()
    end
    
    function update_νonly(val)
        ν = deg2rad(sν.val)
        a = ustrip(u"Rsun", sa.val*u"AU")
        e = se.val
        i = deg2rad(si.val)
        ω = deg2rad(sω.val)
        Ω = deg2rad(sΩ.val)
        
        x, y, z = pos(ν, a, e, i, ω, Ω)
        sec.set_center((x, y))
        sec.set_zorder(z)
        
        νlines(_νlines, ν, ω)
        fig.canvas.draw_idle()
    end
    
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

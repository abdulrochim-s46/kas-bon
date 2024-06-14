// concrete
var fc; // MPa, spesified compressive strength of concrete
var fy_; // MPa, specified yield strength of steel reinforcement
var Ec = 4700.0 * Math.sqrt(fc); // MPa, modulus of elasticity of concrete
var betha_1;
if (fc <= 27.58) { betha_1 = 0.85;}
else { betha_1 = Math.max(0.65, 1.09 - 0.008 * fc);}

// metaldeck
// properties
var Fy; //MPa, specified yield strength of steel deck
var Es; //MPa, modulus of elasticity of steel deck
var As; // mm2/m, area of steel deck per unit slab width (b = 1000 mm)
var Isf; //mm4, moment of inertia of the full (unreduced) steel deck per unit slab width

// dimension
var ln; //mm, length of slab
var h; // mm, over all slab depth
var dd; //mm, depth of deck
var hc; //mm, depth of concrete above steel deck
var Cs; //mm, pitch of deck ribs
var Wr; //mm, average deck rib width
var d; //mm, distance from top of concrete (extreme compression fiber) to centroid of steel deck
var b = 1000.0; //mm, unit width of compression face of composite slab (unit slab width = 1000 mm)

// loads data
var ms; // kg/m2, mass of self slab concrete
var sdl; // kg/m2, superimposed dead loads
var ll; // kg/m2, live loads

// under reinforced or over reinforced
var c_div_d; //compression depth ratio
var c_div_d_balance; //compression depth ratio of balanced condition

//For Cracked Section
var rho; // steel deck ratio
var n; // modular ratio
var ycc_crack; //mm, distance from top of concrete slab to neutral axis of transformed composite section
var ycs_crack; //mm, distance from neutral axis of transformed composite section to centroid of steel deck

//Calculation of the cracked of moment of inertia
var Icr; //mm4, Cracked moment of inertia

//For Uncracked section
var ycc_uncrack; //mm, distance from top of concrete slab to neutral axis of transformed composite section ==> Uncraked Section
var ycs_uncrack; //mm, distance from neutral axis of transformed composite section to centroid of steel deck ==> Uncracked Section

//Calculation of the uncracked of moment of inertia
var I_uncrack; //mm4, Uncracked moment of inertia

//Moment of inertia of the composite section effective
var Id; // mm4

//Flexural Strength
//Under-reinforced Slabs (c/d < c/d balance) => positive moment bending
var phi_s = 0.85;
var My; // Nmm, Yield moment for the composite dec-slab, considering a craked cross section
var Mru;

// Over-reinforced Slabs (c/d >= s/d balance) =>  positive moment bending
var eps_cu = 0.003;
var phi_c = 0.65;
var m;
var c;
var Mro;

//loads combination
var comb1;

// moment ultimate
var Mu; // Nmm/m

//declaration of string safe and sign
var safe;
var sign;
var phi_mn;
var underreinforced_or_overreinforced;
var sign_under_or_over_reinforced;

function steeldeck_design(){
    // concrete
    fc = parseFloat(document.getElementById('fc').value); // MPa, spesified compressive strength of concrete
    fy_ = parseFloat(document.getElementById('fyc').value); // MPa, specified yield strength of steel reinforcement
    Ec = 4700.0 * Math.sqrt(fc); // MPa, modulus of elasticity of concrete
    
    if (fc <= 27.58) { betha_1 = 0.85;}
    else { betha_1 = Math.max(0.65, 1.09 - 0.008 * fc);}

    // metaldeck
    // properties
    Fy = parseFloat(document.getElementById('Fy').value); //MPa, specified yield strength of steel deck
    Es = parseFloat(document.getElementById('Es').value); //MPa, modulus of elasticity of steel deck
    As = parseFloat(document.getElementById('As').value); // mm2/m, area of steel deck per unit slab width (b = 1000 mm)
    Isf = parseFloat(document.getElementById('Isf').value); //mm4, moment of inertia of the full (unreduced) steel deck per unit slab width

    // dimension
    ln = parseFloat(document.getElementById('ln').value); //mm, length of slab
    h = parseFloat(document.getElementById('h').value); // mm, over all slab depth
    dd = parseFloat(document.getElementById('dd').value); //mm, depth of deck
    hc = h - dd; //mm, depth of concrete above steel deck
    Cs = parseFloat(document.getElementById('Cs').value); //mm, pitch of deck ribs
    Wr = parseFloat(document.getElementById('Wr').value); //mm, average deck rib width
    d = h - dd/2.0; //mm, distance from top of concrete (extreme compression fiber) to centroid of steel deck

    // loads data
    ms = 2400.0 * h/1000.0; // kg/m2, mass of self slab concrete
    sdl = parseFloat(document.getElementById('sdl').value); // kg/m2, superimposed dead loads
    ll = parseFloat(document.getElementById('ll').value); // kg/m2, live loads

    // under reinforced or over reinforced
    c_div_d = As * Fy / (fc * d * b * betha_1); //compression depth ratio
    c_div_d_balance = 0.003 * (h - dd)/ ((Fy/ Es + 0.003) * d); //compression depth ratio of balanced condition

    //For Cracked Section
    rho = As/ (b * d); // steel deck ratio
    n = Es/ Ec; // modular ratio
    ycc_crack = Math.min(hc, d * (Math.sqrt(2.0 * rho * n + Math.pow(rho * n, 2.0)) - rho * n)); //mm, distance from top of concrete slab to neutral axis of transformed composite section
    ycs_crack = d - ycc_crack; //mm, distance from neutral axis of transformed composite section to centroid of steel deck

    //Calculation of the cracked of moment of inertia
    Icr = b/(3.0*n) * Math.pow(ycc_crack, 3.0) + As * Math.pow(ycs_crack, 2.0) + Isf; //mm4, Cracked moment of inertia

    //For Uncracked section
    ycc_uncrack = ( 0.5 * b * Math.pow(hc, 2.0) + n * As * d + Wr * dd * (h - 0.5 * dd) * b/ Cs) / ( b * hc + n * As + Wr * dd * b / Cs); //mm, distance from top of concrete slab to neutral axis of transformed composite section ==> Uncraked Section
    ycs_uncrack = d - ycc_uncrack; //mm, distance from neutral axis of transformed composite section to centroid of steel deck ==> Uncracked Section

    //Calculation of the uncracked of moment of inertia
    I_uncrack =  b * Math.pow(hc, 3.0)/ (12.0 * n) + b * hc/ n * Math.pow((ycc_uncrack - 0.5 * hc), 2.0) + Isf + As * Math.pow(ycs_uncrack, 2.0) + Wr * b * dd/ (n * Cs) * ( Math.pow(dd, 2.0)/ 12.0 + Math.pow((h - ycc_uncrack - 0.5 * dd), 2.0) ); //mm4, Uncracked moment of inertia

    //Moment of inertia of the composite section effective
    Id = (Icr + I_uncrack)/ 2.0; // mm4

    //Flexural Strength
    //Under-reinforced Slabs (c/d < c/d balance) => positive moment bending
    My = Fy * Icr / (h - ycc_crack); // Nmm, Yield moment for the composite dec-slab, considering a craked cross section
    Mru = phi_s * My;

    // Over-reinforced Slabs (c/d >= s/d balance) =>  positive moment bending
    m = Es * eps_cu / (fc * betha_1);
    c = d * (Math.sqrt(rho * m + Math.pow(rho * m/2.0, 2.0)) - rho * m/ 2.0);
    Mro = phi_c * fc * b * betha_1 * c * (d - betha_1 * c/ 2.0);

    //loads combination
    comb1 = (1.2 * (ms + sdl) + 1.6 * ll) * 9.81;

    // moment ultimate
    Mu = 1.0/ 14.0 * (comb1) * ln/1000. * ln/1000. * 1000.; // Nmm/m

    if (c_div_d <= c_div_d_balance){
        phi_mn = Mru;
        underreinforced_or_overreinforced = "under-reinforced";
        sign_under_or_over_reinforced = "<";
    }else {
        phi_mn = Mro;
        underreinforced_or_overreinforced = "over-reinforced";
        sign_under_or_over_reinforced = ">";
    }

    if(phi_mn >= Mu){
        safe = "AMAN";
        sign = ">";
    } else{
        safe ="TIDAK AMAN";
        sign = "<";
    }
    // document.getElementById('result').value = convertedValue.toFixed(0) + ' ' + resultUnit;
    document.getElementById('result').value = underreinforced_or_overreinforced +
                                             " ==> phi Mn = " + (phi_mn/1000000).toFixed(2) + " kNm " + sign + " Mu = " + (Mu/1000000).toFixed(2) + " kNm, [" + safe + "]";
    
}


function metalDeckCalculation(){
    var textarea =
    document.getElementById('textarea');
    textarea.value = 
    "LAPORAN PERHITUNGAN\n" +
    "'Berdasarkan American National Standards Institute(ANSI) / Steel Deck Institute(SDI)'\n\n" +
    "dimensi slab:\n" +
    "Ln = " + ln + " mm\n" +
    "h = " + h + " mm\n\n" +

    "material properti beton:\n" +
    "f'c = " + fc + " MPa\n" +
    "fy = " + fy_ + " MPa\n\n" +

    "material properti steel-deck:\n" +
    "Fy = " + Fy + " MPa\n" +
    "Es = " + Es + " MPa\n" +
    "As = " + As + " mm2\n" +
    "Isf = " + Isf + " mm4\n\n" +

    "dimensi steel-deck:\n" +
    "dd = " + dd + " mm\n" +
    "Cs = " + Cs + " mm\n" +
    "Wr = " + Wr + " mm\n\n" +

    "data beban:\n" +
    "beban pelat = " + ms + " kg/m2\n" +
    "beban SDL = " + sdl + " kg/m2\n" +
    "beban LL = " + ll + " kg/m2\n\n" +

    "perhitungan kondisi under-reinforced or over-reinforced?\n" +
    "c/d = " + c_div_d.toFixed(4) + " " + sign_under_or_over_reinforced + " c/d balance = " + c_div_d_balance.toFixed(4) + ", sehingga " + underreinforced_or_overreinforced + " condition\n\n" +

    "hasil kalkulasi:\n" +
    "momen nominal, phi Mn = phi. My\n" +
    "phi Mn = " + (phi_mn/1000000).toFixed(2) + " kN.m/ m' " + sign + " Mu = " + (Mu/1000000).toFixed(2) + " kNm/ m', [" + safe + "]"; 

}


// console.log("ms = " + ms)
// console.log("Ec = " + Ec)
// console.log("hc = " + hc)
// console.log("d = " + d)
// console.log("betha_1 = " + betha_1)
// console.log("c/d = " + c_div_d)
// console.log("c/d balance= " + c_div_d_balance)
// console.log("ycc_crack = " + ycc_crack)
// console.log("ycs_crack = " + ycs_crack)
// console.log("Icr = " + Icr)

// console.log("ycc_uncrack = " + ycc_uncrack)
// console.log("ycs_uncrack = " + ycs_uncrack)
// console.log("I_uncrack = " + I_uncrack)
// console.log("I_eff = " + Id)
// console.log("Mru = " + Mru)
// console.log("c/d = " + c/d)
// console.log("Mro = " + Mro)

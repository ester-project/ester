#include "ester-config.h"
#include "utils.h"
#include "matrix.h"
#include "star.h"
#include <string.h>
#include <cstdlib>
#ifdef USE_HDF5
#include <H5Cpp.h>
#endif
#include "matplotlib.h"

matrix star2d::solve_phi() {
    symbolic S;
    S.set_map(map);

    sym symPhi = S.regvar("Phi");
    matrix phi, rho;

    phi = zeros(nr, nth);
    rho = this->rho;

    sym eq = lap(symPhi);
    solver op;
    op.init(ndomains, 1, "full");
    op.regvar("Phi");
    op.set_nr(map.npts);

    S.set_value("Phi", phi);

    op.reset();
    op.set_nr(map.npts);

    matrix rhs = rho*pi_c;

    eq.add(&op, "Phi", "Phi");

    int j0 = 0;
    for (int n=0; n<ndomains; n++) {
        if (n == 0) {
            // BC at center
            op.bc_bot2_add_l(0, "Phi", "Phi", ones(1, nth), map.D.block(0).row(0));
            rhs(0) = 0.0;
        }
        else {
            op.bc_bot2_add_l(n, "Phi", "Phi", ones(1,1), D.block(n).row(0));
            op.bc_bot1_add_l(n, "Phi", "Phi", -ones(1,1), D.block(n-1).row(-1));
            rhs(j0) = -(D, phi)(j0) + (D, phi)(j0-1);
        }

        if (n == ndomains-1) {
            // BC at surface
            op.bc_top1_add_l(n, "Phi", "Phi", ones(1, nth), map.D.block(n).row(-1));
            op.bc_top1_add_d(n, "Phi", "Phi", ones(1, nth));
            rhs(-1) = -(map.D, phi)(-1) - phi(-1);
        }
        else {
            op.bc_top1_add_d(n,"Phi", "Phi", ones(1,1));
            op.bc_top2_add_d(n,"Phi", "Phi", -ones(1,1));
            rhs(j0+map.gl.npts[n]-1) = -phi(j0+map.gl.npts[n]-1)+phi(j0+map.gl.npts[n]);
        }
        j0 += map.gl.npts[n];
    }

    op.set_rhs("Phi", rhs);

    op.solve();

    phi = op.get_var("Phi");

    return phi;
}

star2d::star2d() : nr(map.gl.N), nth(map.leg.npts), nex(map.ex.gl.N),
    ndomains(map.gl.ndomains), r(map.r), z(map.z), th(map.th), Dt(map.Dt),
    Dt2(map.Dt2), zex(map.ex.z), Dex(map.ex.D), rex(map.ex.r),
    D(map.D) {

    config.newton_dmax=0.5;
    config.verbose=0;
    version.name = std::string(VERSION);

    const char *minor = strstr(VERSION, ".") + 1;
    const char *rev = strstr(minor, ".") + 1;
    version.major=atoi(VERSION);
    version.minor=atoi(minor);
    version.rev=atoi(rev);
    version.svn=0;
    stratified_comp = 0;
    config.dump_iter = 0;
}

star2d::~star2d() {
}

star2d::star2d(const star2d &A) : nr(map.gl.N), nth(map.leg.npts),
    nex(map.ex.gl.N), ndomains(map.gl.ndomains), r(map.r), z(map.z), th(map.th),
    Dt(map.Dt), Dt2(map.Dt2), zex(map.ex.z), Dex(map.ex.D), rex(map.ex.r),
    D(map.D) {

        copy(A);
    }

star2d &star2d::operator=(const star2d &A) {

    copy(A);

    return *this;
}

void star2d::copy(const star2d &A) {

    phi=A.phi;
    p=A.p;
    T=A.T;
    rho=A.rho;
    comp=A.comp;

    opa=A.opa;
    eos=A.eos;
    nuc=A.nuc;
    atm=A.atm;
    config=A.config;
    units=A.units;

    R=A.R;M=A.M;
    Tc=A.Tc;pc=A.pc;rhoc=A.rhoc;
    X0=A.X0;Y0=A.Y0;Z0=A.Z0;
    surff=A.surff;
    conv=A.conv;
    Xc=A.Xc;

    map=A.map;

    ps=A.ps;Ts=A.Ts;
    m=A.m;pi_c=A.pi_c;Lambda=A.Lambda;

    phiex=A.phiex;
    w=A.w;
    vr=A.vr;vt=A.vt;G=A.G;

    Omega=A.Omega;Omega_bk=A.Omega_bk;
    Ekman=A.Ekman;

    core_convec=A.core_convec;
    env_convec=A.env_convec;
    min_core_size=A.min_core_size;

    domain_type=A.domain_type;

}

#ifdef USE_HDF5
template<typename T>
void write_attr(H5::Group grp, const char *name, H5::DataType type,
        const T ptr, int n = 1) {
    H5::DataSpace dataspace;
    if (n == 1) {
        dataspace = H5::DataSpace(H5S_SCALAR);
    }
    else {
        hsize_t dims[1];
        dims[0] = n;
        dataspace = H5::DataSpace(1, dims);
    }
    H5::Attribute attr = grp.createAttribute(name, type, dataspace);
    attr.write(type, ptr);
}

void write_field(H5::Group grp, const char *name, const matrix &field) {
    hsize_t dims[2];
    dims[0] = field.ncols();
    dims[1] = field.nrows();
    H5::DataSpace dataspace(2, dims);
    H5::DataSet dataset = grp.createDataSet(name, H5::PredType::IEEE_F64LE,
            dataspace);
    dataset.write(field.data(), H5::PredType::IEEE_F64LE);
}
#endif

void star2d::hdf5_write(const char *filename) const {
#ifdef USE_HDF5
#ifndef DEBUG
    H5::Exception::dontPrint();
#endif

    H5::H5File file(filename, H5F_ACC_TRUNC);
    H5::Group star = file.createGroup("/star");

    H5::PredType integer = H5::PredType::STD_I32LE;
    H5::PredType real    = H5::PredType::IEEE_F64LE;

    write_attr(star, "ndomains",    integer,    &map.ndomains);
    write_attr(star, "nr",          integer,    &nr);
    write_attr(star, "npts",        integer,    map.gl.npts, map.ndomains);
    write_attr(star, "nth",         integer,    &map.leg.npts);
    write_attr(star, "nex",         integer,    &map.ex.gl.npts[0]);
    write_attr(star, "conv",        integer,    &conv);
    write_attr(star, "domain_type", integer,    domain_type.data(), map.ndomains);
    write_attr(star, "core_convec", integer,    &core_convec);
    write_attr(star, "env_convec",  integer,    &env_convec);
    write_attr(star, "stratified_comp", integer, &stratified_comp);
    write_attr(star, "xif",         real,       map.gl.xif, map.ndomains+1);
    write_attr(star, "M",           real,       &M);
    write_attr(star, "R",           real,       &R);
    write_attr(star, "X0",          real,       &X0);
    write_attr(star, "Z0",          real,       &Z0);
    write_attr(star, "Xc",          real,       &Xc);
    write_attr(star, "surff",       real,       &surff);
    write_attr(star, "Tc",          real,       &Tc);
    write_attr(star, "pc",          real,       &pc);
    write_attr(star, "Omega",       real, &Omega);
    write_attr(star, "Omega_bk",    real, &Omega_bk);
    write_attr(star, "Ekman",       real, &Ekman);
    write_attr(star, "min_core_size", real, &min_core_size);
    double ek = this->virial_L()/2;
    write_attr(star, "Ek",          real, &ek);

    H5::StrType strtype;
    strtype = H5::StrType(H5::PredType::C_S1, strlen(version.name.c_str())+1);
    write_attr(star, "version", strtype, H5std_string(version.name));

    strtype = H5::StrType(H5::PredType::C_S1, strlen(opa.name)+1);
    write_attr(star, "opa.name", strtype, H5std_string(opa.name));

    strtype = H5::StrType(H5::PredType::C_S1, strlen(eos.name)+1);
    write_attr(star, "eos.name", strtype, eos.name);

    strtype = H5::StrType(H5::PredType::C_S1, strlen(nuc.name)+1);
    write_attr(star, "nuc.name", strtype, nuc.name);

    strtype = H5::StrType(H5::PredType::C_S1, strlen(atm.name)+1);
    write_attr(star, "atm.name", strtype, atm.name);

    matrix_map fields;
    fields["r"] = r;
    fields["z"] = z;
    fields["th"] = th;
    fields["rho"] = rho;
    fields["phi"] = phi;
    fields["phiex"] = phiex;
    fields["R"] = map.R;
    fields["p"] = p;
    fields["T"] = T;
    fields["G"] = G;
    fields["w"] = w;
    fields["X"] = comp.X();
    fields["Y"] = comp.Y();
    fields["Z"] = comp.Z();
    fields["N2"] = N2();
    fields["nuc.eps"] = nuc.eps;

    for (matrix_map::iterator it=fields.begin(); it!=fields.end(); ++it) {
        write_field(star, it->first.c_str(), it->second);
    }
#endif
}

void star2d::write(const char *output_file, char mode) const {
    OUTFILE fp;

    if (isHDF5Name(output_file)) {
        hdf5_write(output_file);
        return;
    }

    fp.open(output_file,mode);
    write_tag(&fp);

    fp.write("ndomains",&ndomains);
    fp.write("npts",map.gl.npts,ndomains);
    fp.write("xif",map.gl.xif,ndomains+1);
    fp.write("nth",&map.leg.npts,1);
    fp.write("nex",map.ex.gl.npts,1);
    fp.write("M",&M);
    fp.write("R",&R);
    fp.write("X0",&X0);
    fp.write("Z0",&Z0);
    fp.write("Xc",&Xc);
    fp.write("conv",&conv);
    fp.write("domain_type",&domain_type[0],ndomains);
    fp.write("surff",&surff);
    fp.write("Tc",&Tc);
    fp.write("pc",&pc);
    fp.write("opa.name",opa.name,strlen(opa.name)+1);
    fp.write("eos.name",eos.name,strlen(eos.name)+1);
    fp.write("nuc.name",nuc.name,strlen(nuc.name)+1);
    fp.write("atm.name",atm.name,strlen(atm.name)+1);
    fp.write("Omega",&Omega);
    fp.write("Omega_bk",&Omega_bk);
    fp.write("Ekman",&Ekman);
    fp.write("core_convec",&core_convec);
    fp.write("env_convec",&env_convec);
    fp.write("min_core_size",&min_core_size);
    fp.write("stratified_comp",&stratified_comp);
    fp.write("version.major",&version.major);
    fp.write("version.minor",&version.minor);
    fp.write("version.rev",&version.rev);
    fp.write("version.svn",&version.svn);

    fp.write("phi",&phi);
    fp.write("p",&p);
    fp.write("T",&T);
    fp.write("phiex",&phiex);
    fp.write("map.R",&map.R);
    fp.write("w",&w);
    fp.write("G",&G);
    fp.write("comp",(matrix_map *)&comp);

    fp.close();

}

#ifdef USE_HDF5
template<typename T>
int read_attr(H5::Group grp, const char *name, T mem) {
    try {
        H5::Attribute attr = grp.openAttribute(name);
        H5::DataType type = H5::DataType(attr.getDataType());
        attr.read(type, mem);
    }
    catch (H5::Exception e) {
        return 1;
    }
    return 0;
}

int read_field(H5::Group grp, const char *name, matrix &field) {
    try {
        H5::DataSet set = grp.openDataSet(name);
        H5::DataSpace filespace = set.getSpace();
        hsize_t dims[2];
        if (filespace.getSimpleExtentDims(dims, NULL) != 2) {
            return 1;
        }
        H5::DataSpace mspace(2, dims);
        field.dim(dims[1], dims[0]);
        set.read(field.data(), H5::PredType::IEEE_F64LE, mspace, filespace);
    }
    catch (H5::Exception) {
        return 1;
    }
    return 0;
}
#endif

int star2d::hdf5_read(const char *input_file, int dim) {
#ifdef USE_HDF5
#ifndef DEBUG
    H5::Exception::dontPrint();
#endif

    H5::H5File file;
    H5std_string buf;
    try {
        file = H5::H5File(input_file, H5F_ACC_RDONLY);
    }
    catch (H5::FileIException e) {
        ester_err("Could not open file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    H5::Group star;
    try {
        star = file.openGroup("/star");
    }
    catch (H5::Exception) {
        ester_err("Could not open '/star' in `%s'", input_file);
        exit(EXIT_FAILURE);
    }

    int ndoms;
    if (read_attr(star, "nth", &map.leg.npts)) {
        ester_err("could not read 'nth' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    if ((map.leg.npts == 1 && dim == 2) || (map.leg.npts > 1 && dim == 1)) {
        return 1;
    }
    if (read_attr<H5std_string&>(star, "version", buf)) {
        ester_warn("Could not read 'version' from file `%s'", input_file);
        version.name = "unknown";
    }
    version.name = buf;
    if (read_attr(star, "ndomains", &ndoms)) {
        ester_err("Could not read 'ndomains' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    map.gl.set_ndomains(ndoms);
    if (read_attr(star, "npts", map.gl.npts)) {
        ester_err("Could not read 'npts' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    if (read_attr(star, "xif", &map.gl.xif[0])) {
        ester_err("Could not read 'xif' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    if (read_attr(star, "nex", map.ex.gl.npts)) {
        ester_err("Could not read 'nex' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    if (read_attr(star, "M", &M)) {
        ester_warn("Could not read 'M' from file `%s'", input_file);
    }
    if (read_attr(star, "R", &R)) {
        ester_warn("Could not read 'R' from file `%s'", input_file);
    }
    if (read_attr(star, "X0", &X0)) {
        ester_warn("Could not read 'X0' from file `%s'", input_file);
    }
    if (read_attr(star, "Z0", &Z0)) {
        ester_warn("Could not read 'Z0' from file `%s'", input_file);
    }
    if (read_attr(star, "Xc", &Xc)) {
        ester_warn("Could not read 'Xc' from file `%s'", input_file);
    }
    if (read_attr(star, "conv", &conv)) {
        ester_err("Could not read 'conv' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }

    domain_type.resize(ndomains);
    if (read_attr(star, "domain_type", &domain_type[0])) {
        for (int n=0; n<ndomains; n++) {
            if (n < conv) domain_type[n] = CORE;
            else domain_type[n] = RADIATIVE;
        }
    }
    if (read_attr(star, "surff", &surff)) {
        ester_warn("Could not read 'surff' from file `%s'", input_file);
    }
    if (read_attr(star, "Tc", &Tc)) {
        ester_warn("Could not read 'Tc' from file `%s'", input_file);
    }
    if (read_attr(star, "pc", &pc)) {
        ester_warn("Could not read 'pc' from file `%s'", input_file);
    }

    if (read_attr<H5std_string&>(star, "opa.name", buf)) {
        ester_warn("Could not read 'opa.name' from file `%s'", input_file);
        buf = H5std_string("opal");
    }
    strncpy(opa.name, buf.c_str(), 15);

    if (read_attr<H5std_string&>(star, "eos.name", buf)) {
        ester_warn("Could not read 'eos.name' from file `%s'", input_file);
        buf = H5std_string("opal");
    }
    strncpy(eos.name, buf.c_str(), 15);

    if (read_attr<H5std_string&>(star, "nuc.name", buf)) {
        ester_warn("Could not read 'nuc.name' from file `%s'", input_file);
        buf = H5std_string("simple");
    }
    strncpy(nuc.name, buf.c_str(), 15);

    if (read_attr<H5std_string&>(star, "atm.name", buf)) {
        ester_warn("Could not read 'atm.name' from file `%s'", input_file);
        buf = H5std_string("simple");
    }
    strncpy(atm.name, buf.c_str(), 15);

    if (read_attr(star, "Omega", &Omega)) {
        ester_warn("Could not read 'Omega' from file `%s'", input_file);
        Omega = .0;
    }
    if (read_attr(star, "Omega_bk", &Omega_bk)) {
        ester_warn("Could not read 'Omega_bk' from file `%s'", input_file);
        Omega_bk = .0;
    }
    if (read_attr(star, "Ekman", &Ekman)) {
        ester_warn("Could not read 'Ekman' from file `%s'", input_file);
        Ekman = .0;
    }
    if (read_attr(star, "core_convec", &core_convec)) {
        ester_warn("Could not read 'core_convec' from file `%s'", input_file);
        core_convec = 1;
    }
    if (read_attr(star, "env_convec", &env_convec)) {
        ester_warn("Could not read 'env_convec' from file `%s'", input_file);
        env_convec = 0;
    }
    if (read_attr(star, "min_core_size", &min_core_size)) {
        ester_warn("Could not read 'min_core_size' from file `%s'", input_file);
        min_core_size = 0.03;
    }
    if (read_attr(star, "stratified_comp", &stratified_comp)) {
        stratified_comp = 0;
    }

    map.init();

    if (read_field(star, "phi", phi)) {
        ester_err("Could not read field 'phi' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    if (read_field(star, "p", p)) {
        ester_err("Could not read field 'p' from file `%s'", input_file);
        exit(EXIT_FAILURE);
    }
    if (read_field(star, "T", T)) {
        ester_err("Could not read field 'T' from file `%s'", input_file);
        T = zeros(nr, nth);
        exit(EXIT_FAILURE);
    }
    if (read_field(star, "phiex", phiex)) {
        ester_warn("Could not read field 'phiex' from file `%s'", input_file);
        phiex = zeros(nr, nth);
    }
    if (read_field(star, "R", map.R)) {
        ester_err("Could not read field 'R' from file `%s'", input_file);
        map.R = zeros(nr, nth);
        exit(EXIT_FAILURE);
    }
    if (map.R.nrows() < map.ndomains+1)
        map.R = zeros(1,nth).concatenate(map.R);
    map.remap();
    if (read_field(star, "w", w)) {
        ester_warn("Could not read field 'w' from file `%s'", input_file);
        w = zeros(nr, nth);
    }
    if (read_field(star, "G", G)) {
        ester_warn("Could not read field 'G' from file `%s'", input_file);
        G = zeros(nr, nth);
    }
    if (read_field(star, "X", comp["H"])) {
        ester_warn("Could not read field 'X' from file `%s'", input_file);
        comp["H"] = zeros(nr, nth);
    }

    fill();

    matrix sol_phi = solve_phi();
    phi = sol_phi;

    return 0;
#else
    ester_err("Could not read from hdf5 file: HDF5 is not enabled");
#endif
}

int star2d::read(const char *input_file, int dim) {
    char tag[32];
    int ndom;
    INFILE fp;

    // if input file ends with '.hdf5': read in hdf5 format
    if (isHDF5Name(input_file)) {
        return hdf5_read(input_file, dim);
    }

    if(!fp.open(input_file,'b')) {
        if(!fp.open(input_file,'t')) {
            return read_old(input_file);
        }
    }

    tag[0]='\0';
    if(fp.len("tag")<=16) fp.read("tag",tag);

    tag[16]='\0';
    if(!check_tag(tag)) {
        fp.close();
        return 1;
    }

    if(fp.read("version.major",&version.major)) {
        char ver[32];
        if(fp.read("version",ver)) {
            version.major=1;
            version.minor=0;
            version.rev=70;
            if(!strcmp(ver,"1.0.0")) version.svn=0;
            else version.svn=1;
        } else {
            version.major=0;
            version.minor=0;
            version.rev=0;
            version.svn=1;
        }
    } else {
        fp.read("version.minor",&version.minor);
        fp.read("version.rev",&version.rev);
        fp.read("version.svn",&version.svn);
    }
    char *buf = NULL;
    if (version.svn) {
        if (asprintf(&buf, "%d.%d rev %d",
                    version.major,
                    version.minor,
                    version.rev) == -1) {
            ester_err("out of memory");
        }
    }
    else {
        if (asprintf(&buf, "%d.%d.%d",
                    version.major,
                    version.minor,
                    version.rev) == -1) {
            ester_err("out of memory");
        }
    }
    version.name = std::string(buf);
    free(buf);
    fp.read("ndomains",&ndom);
    map.gl.set_ndomains(ndom);
    fp.read("npts",map.gl.npts);
    fp.read("xif",map.gl.xif);
    if(fp.read("nth",&map.leg.npts)) map.leg.npts=1;
    fp.read("nex",map.ex.gl.npts);
    fp.read("M",&M);
    fp.read("R",&R);
    if(fp.read("X0",&X0)) fp.read("X",&X0);
    if(fp.read("Z0",&Z0)) fp.read("Z",&Z0);
    fp.read("Xc",&Xc);
    fp.read("conv",&conv);
    domain_type.resize(ndomains);
    if(fp.read("domain_type",&domain_type[0])) {
        for(int n=0;n<ndomains;n++) {
            if(n<conv) domain_type[n]=CORE;
            else domain_type[n]=RADIATIVE;
        }
    }
    fp.read("surff",&surff);
    fp.read("Tc",&Tc);
    fp.read("pc",&pc);
    fp.read("opa.name",opa.name);
    fp.read("eos.name",eos.name);
    fp.read("nuc.name",nuc.name);
    if(fp.read("atm.name",atm.name)) strcpy(atm.name,"simple");
    if(fp.read("Omega",&Omega)) Omega=0;
    if(fp.read("Omega_bk",&Omega_bk)) Omega_bk=0;
    if(fp.read("Ekman",&Ekman)) Ekman=0;
    if(fp.read("core_convec",&core_convec)) core_convec=1;
    if(fp.read("env_convec",&env_convec)) env_convec=0;
    if(fp.read("min_core_size",&min_core_size)) min_core_size=0.03;
    if(fp.read("stratified_comp",&stratified_comp)) stratified_comp = 0;

    map.init();

    fp.read("phi",&phi);
    fp.read("p",&p);
    fp.read("T",&T);
    if(fp.read("phiex",&phiex)) phiex=zeros(nex,nth);
    fp.read("map.R",&map.R);
    if(map.R.nrows()<map.ndomains+1) {
        map.R=zeros(1,nth).concatenate(map.R);
    }
    map.remap();
    if(fp.read("w",&w)) {
        w=zeros(nr,nth);
        ester_warn("cannot read w,set to 0\n");
    }
    if(fp.read("G",&G)) G=zeros(nr,nth);
    if(fp.read("comp",(matrix_map *)&comp)) init_comp();

    fp.close();
    fill();

    return 0;
}

void star2d::write_tag(OUTFILE *fp) const {
    char tag[7]="star2d";

    fp->write("tag",tag,7);

}

bool star2d::check_tag(const char *tag) const {
    if(strcmp(tag,"star2d")) return false;
    return true;

}

int star2d::read_old(const char *input_file){
    FILE *fp;
    char tag[7],mode,*c;
    int ndom;

    if(!(fp=fopen(input_file,"rb"))) {
        return 0;
    }
    if (fread(tag,1,7,fp) != 7)
        ester_warn("read failed");
    tag[6]='\0';
    if(strcmp(tag,"star2d")) {
        return 0;
    }
    if (fread(&mode,1,1,fp) != 1)
        ester_warn("read failed");
    fclose(fp);

    if(mode=='b')
        fp=fopen(input_file,"rb");
    else
        fp=fopen(input_file,"rt");

    if (fread(tag,1,7,fp) != 7)
        ester_warn("read failed");
    if (fread(&mode,1,1,fp) != 1)
        ester_warn("read failed");

    if(mode=='b') {
        if (fread(&ndom,sizeof(int),1,fp) != 1) ester_warn("read failed");
        map.gl.set_ndomains(ndom);
        if (fread(map.gl.npts,sizeof(int),ndom,fp) != (size_t) ndom)
            ester_warn("read failed");
        if (fread(map.gl.xif,sizeof(double),ndom+1,fp) != (size_t) ndom+1)
            ester_warn("read failed");
        if (fread(map.ex.gl.npts,sizeof(int),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&map.leg.npts,sizeof(int),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&M,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&R,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&X0,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&Z0,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&Xc,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&conv,sizeof(int),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&surff,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&Omega,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&Omega_bk,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&Tc,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        if (fread(&pc,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
        c=opa.name;
        *c=fgetc(fp);
        while(*(c++)!='\0') *c=fgetc(fp);
        c=eos.name;
        *c=fgetc(fp);
        while(*(c++)!='\0') *c=fgetc(fp);
        c=nuc.name;
        *c=fgetc(fp);
        while(*(c++)!='\0') *c=fgetc(fp);
        c=atm.name;
        *c=fgetc(fp);
        while(*(c++)!='\0') *c=fgetc(fp);
    } else {
        int i;
        if (fscanf(fp,"\n%d ",&ndom) != 1)
            ester_warn("read failed");
        map.gl.set_ndomains(ndom);
        for(i=0;i<ndom;i++)
            if (fscanf(fp,"%d ",(map.gl.npts+i)) != 1)
                ester_warn("read failed");
        for(i=0;i<ndom+1;i++)
            if (fscanf(fp,"%le ",(map.gl.xif+i)) != 1)
                ester_warn("read failed");
        if (fscanf(fp,"%d %d",map.ex.gl.npts,&map.leg.npts) != 2)
            ester_warn("read failed");
        if (fscanf(fp,"\n%le %le %le %le\n",&M,&R,&X0,&Z0) != 4)
            ester_warn("read failed");
        if (fscanf(fp,"%le %d %le\n",&Xc,&conv,&surff) != 3)
            ester_warn("read failed");
        if (fscanf(fp,"%le %le\n",&Omega,&Omega_bk) != 2)
            ester_warn("read failed");
        if (fscanf(fp,"%le %le\n",&Tc,&pc) != 2)
            ester_warn("read failed");
        if (fscanf(fp,"%s\n",opa.name) != 1)
            ester_warn("read failed");
        if (fscanf(fp,"%s\n",eos.name) != 1)
            ester_warn("read failed");
        if (fscanf(fp,"%s\n",nuc.name) != 1)
            ester_warn("read failed");
        if (fscanf(fp,"%s\n",atm.name) != 1)
            ester_warn("read failed");
    }
    map.init();
    map.R.read(ndomains,nth,fp,mode);
    map.R=zeros(1,nth).concatenate(map.R);
    map.remap();
    phi.read(nr,nth,fp,mode);
    phiex.read(nex,nth,fp,mode);
    p.read(nr,nth,fp,mode);
    T.read(nr,nth,fp,mode);
    w.read(nr,nth,fp,mode);
    G.read(nr,nth,fp,mode);
    matrix psi;
    psi.read(nr,nth,fp,mode);
    if(mode=='b') {
        if (fread(&Ekman,sizeof(double),1,fp) != 1)
            ester_warn("read failed");
    } else {
        if (fscanf(fp,"%le\n",&Ekman) != 1)
            ester_warn("read failed");
    }
    core_convec=1;
    env_convec=1;
    min_core_size=0.03;
    version.major=0;
    version.minor=0;
    version.rev=0;
    version.svn=1;
    domain_type.resize(ndomains);
    for(int n=0;n<ndomains;n++) {
        if(n<conv) domain_type[n]=CORE;
        else domain_type[n]=RADIATIVE;
    }
    init_comp();
    fclose(fp);
    fill();
    return 0;
}

int star2d::init(const char *input_file,const char *param_file,int argc,char *argv[]) {
    cmdline_parser cmd;
    file_parser fp;
    char *arg,*val,default_params[256];
    mapping map0;
    int i,k,change_grid=0,nt=-1,next=-1;
    star1d in1d;
    diff_leg leg_new;
    matrix Tr,m0;

    sprintf(default_params,"%s/ester/1d_default.par", ESTER_DATADIR);

    if(input_file != NULL) {
        if (read(input_file)) {
            if(!in1d.read(input_file)) {
                if(*param_file) {
                    if(fp.open(param_file)) {
                        while((k=fp.get(arg,val))) {
                            if(!strcmp(arg,"nth")&&val) nt=atoi(val);
                            if(!strcmp(arg,"nex")&&val) next=atoi(val);
                        }
                        fp.close();
                    }
                }
                cmd.open(argc,argv);
                while(cmd.get(arg,val)) {
                    if(!strcmp(arg,"nth")&&val) nt=atoi(val);
                    if(!strcmp(arg,"nex")&&val) next=atoi(val);
                }
                cmd.close();
                init1d(in1d, nt, next);
            } else {
                ester_err("Error reading input file: %s", input_file);
            }
        }
        map0=map;
    } else {
        if(!fp.open(default_params)) {
            ester_err("Can't open default parameters file %s\n", default_params);
        }
        else {
            while((k=fp.get(arg,val))) {
                if((i=check_arg(arg,val,&change_grid))) {
                    ester_err("Syntax error in parameters file %s, line %d\n",default_params,k);
                    if(i==2) {
                        ester_err("Error: Argument to '%s' missing\n",arg);
                        return 1;
                    }
                    if(i==1) {
                        ester_err("Unknown parameter %s\n",arg);
                        return 1;
                    }
                }
            }
            fp.close();
        }
        change_grid=0;
    }

    if(*param_file) {
        if(!fp.open(param_file)) {
            ester_err("Can't open parameters file %s\n",param_file);
            return 1;
        }
        else {
            while((k=fp.get(arg,val))) {
                if((i=check_arg(arg,val,&change_grid))) {
                    ester_err("Syntax error in parameters file %s, line %d\n",
                            param_file, k);
                    if(i==2) {
                        ester_err("Error: Argument to '%s' missing\n",arg);
                        return 1;
                    }
                    if(i==1) {
                        ester_err("Unknown parameter %s\n",arg);
                        return 1;
                    }
                }
            }
            fp.close();
        }
    }

    cmd.open(argc,argv);
    while(int err_code=cmd.get(arg,val)) {
        if(err_code==-1) exit(1);
        err_code=check_arg(arg,val,&change_grid);
        if(err_code==2) {
            ester_err("Error: Argument to '%s' missing\n",arg);
            return 1;
        }
        if(err_code==1) {
            ester_err("Unknown parameter '%s'\n",arg);
            return 1;
        }
        cmd.ack(arg,val);
    }
    cmd.close();

    if((change_grid&1)&&!(change_grid&2)) {
        ester_err("Must specify number of points per domain (npts)\n");
        exit(1);
    }

    if (*input_file) {
        if(change_grid) {
            mapping map_new;
            map_new=map;
            map=map0;
            remap(map_new.ndomains,map_new.gl.npts,map_new.nt,map_new.nex);
        }
        if(version.rev<=71) { // Force remapping for old files
            int npts[ndomains+1];
            for(int n=0;n<ndomains;n++) npts[n]=map.npts[n];
            npts[ndomains]=npts[ndomains-1];
            remap(ndomains+1,npts,map.nt,map.nex);
            remap(ndomains-1,map.npts,map.nt,map.nex);
        }
    } else {
        map.init();
        T=1-0.5*r*r;
        p=T;
        phi=-T;
        phiex=zeros(nex,nth);
        w=zeros(nr,nth);
        G=zeros(nr,nth);
        conv=0;
        domain_type.resize(ndomains);
        for(int n=0;n<ndomains;n++) domain_type[n]=RADIATIVE;
    }
    init_comp();
    fill();
    return 0;
}

// used when computing a 2D model from a 1D model
void star2d::init1d(const star1d &A,int npts_th,int npts_ex) {
    matrix thd;
    char *arg,*val,default_params[256];
    int k;
    file_parser fp;

    sprintf(default_params,"%s/config/2d_default.par",ESTER_DATADIR);

    *this=A;

    phiex=phi(nr-1)/map.ex.r;
    Omega_bk=0;
    Omega=0;

    if(fp.open(default_params)) {
        while((k=fp.get(arg,val))) {
            if(!strcmp(arg,"nth")&&val&&npts_th==-1) npts_th=atoi(val);
            if(!strcmp(arg,"nex")&&val&&npts_ex==-1) npts_ex=atoi(val);
            if(!strcmp(arg,"Omega_bk")&&val) Omega_bk=atof(val);
            if(!strcmp(arg,"Ekman")&&val) Ekman=atof(val);
        }
        fp.close();
    }

    if(npts_th==-1) npts_th=8;
    if(npts_ex==-1) npts_ex=8;


    remap(ndomains,map.gl.npts,npts_th,npts_ex);

    fill();

}

void star2d::interp(remapper *red) {

    p=red->interp(p);

    phi=red->interp(phi);
    T=red->interp(T);
    w=red->interp(w);
    G=red->interp(G,11);
    comp=red->interp(comp);
    phiex=red->interp_ex(phiex);
    fill(); // recompute the micro-physics variables

}

extern bool dump_jac;
int star2d::check_arg(char *arg,char *val,int *change_grid) {
    int err=0,i;
    char *tok;

    if(!strcmp(arg,"ndomains")) {
        if(val==NULL) return 2;
        map.gl.set_ndomains(atoi(val));
        *change_grid=*change_grid|1;
        if(*change_grid&2) {
            ester_err("ndomains must be specified before npts\n");
            exit(1);
        }
    }
    else if(!strcmp(arg,"npts")) {
        if(val==NULL) return 2;
        tok=strtok(val,",");
        i=0;
        while(tok!=NULL) {
            *(map.gl.npts+i)=atoi(tok);
            tok=strtok(NULL,",");
            i++;
        }
        if(i==1) {
            for(i=1;i<map.gl.ndomains;i++) {
                *(map.gl.npts+i)=*map.gl.npts;
            }
        }
        *change_grid=*change_grid|2;
    }
    else if(!strcmp(arg,"nth")) {
        if(val==NULL) return 2;
        map.leg.npts=atoi(val);
        *change_grid=*change_grid|8;
    }
    else if(!strcmp(arg,"nex")) {
        if(val==NULL) return 2;
        map.ex.gl.set_npts(atoi(val));
        *change_grid=*change_grid|4;
    }
    else if(!strcmp(arg,"M")) {
        if(val==NULL) return 2;
        M=atof(val)*M_SUN;
    }
    else if(!strcmp(arg,"X")) {
        if(val==NULL) return 2;
        X0=atof(val);
    }
    else if(!strcmp(arg,"Z")) {
        if(val==NULL) return 2;
        Z0=atof(val);
    }
    else if(!strcmp(arg,"Xc")) {
        if(val==NULL) return 2;
        Xc=atof(val);
    }
    else if(!strcmp(arg,"conv")) {
        ester_err("Param. conv is no longer modifiable. Disable core convection with core_convec 0.\n");
        return 1;
    }
    else if(!strcmp(arg,"surff")) {
        if(val==NULL) return 2;
        surff=atof(val);
    }
    else if(!strcmp(arg,"Omega_bk")) {
        if(val==NULL) return 2;
        Omega_bk=atof(val);
    }
    else if(!strcmp(arg,"Tc")) {
        if(val==NULL) return 2;
        Tc=atof(val);
    }
    else if(!strcmp(arg,"pc")) {
        if(val==NULL) return 2;
        pc=atof(val);
    }
    else if(!strcmp(arg,"opa")) {
        if(val==NULL) return 2;
        strcpy(opa.name,val);
    }
    else if(!strcmp(arg,"eos")) {
        if(val==NULL) return 2;
        strcpy(eos.name,val);
    }
    else if(!strcmp(arg,"nuc")) {
        if(val==NULL) return 2;
        strcpy(nuc.name,val);
    }
    else if(!strcmp(arg,"atm")) {
        if(val==NULL) return 2;
        strcpy(atm.name,val);
    }
    else if(!strcmp(arg,"Ekman")) {
        if(val==NULL) return 2;
        Ekman=atof(val);
    }
    else if(!strcmp(arg,"core_convec")) {
        if(val==NULL) return 2;
        core_convec=atoi(val);
    }
    else if(!strcmp(arg,"env_convec")) {
        if(val==NULL) return 2;
        env_convec=atoi(val);
    }
    else if(!strcmp(arg,"min_core_size")) {
        if(val==NULL) return 2;
        min_core_size=atof(val);
    }
    else if(!strcmp(arg,"stratified_comp")) {
        if(val==NULL) return 2;
        stratified_comp = atoi(val);
    }
    else if(!strcmp(arg, "dump_jac")) {
        dump_jac = true;
    }
    else if(!strcmp(arg,"dump_iter")) {
        config.dump_iter = 1;
#ifndef USE_HDF5
        ester_warn("ESTER is not built with HDF5 support:\n%s",
                "the -dump_iter option will be inefective");
#endif
    }
    else err=1;

    return err;
}

// dump_info used by 'ester info'
void star2d::dump_info() {
    printf("ESTER 2d model file");
    printf(" (Version %s)\n", version.name.c_str());
    // printf("\n2d ESTER model file (Version %d.%d rev %d",
    //         version.major,
    //         version.minor,
    //         version.rev);
    // if(version.svn) printf(".svn");
    // printf(")\n\n");

    printf("General parameters:\n\n");
    printf("\tMass = %.5f Msun (%e g)\n",M/M_SUN,M);
    printf("\tRadius (p) = %.5f Rsun (%e cm)\n",R/R_SUN,R);
    double re=map.leg.eval_00(r.row(-1),PI/2)(0);
    printf("\tRadius (e) = %.5f Rsun (%e cm)\n",R/R_SUN*re,R*re);
    printf("\tFlatness = %.3f\n",1.-1./re);
    printf("\tLuminosity = %.4f Lsun (%e erg/s)\n",luminosity()/L_SUN,luminosity());
    printf("\tTeff (p) = %.2f\n",map.leg.eval_00(Teff(),0)(0));
    printf("\tTeff (e) = %.2f\n",map.leg.eval_00(Teff(),PI/2)(0));
    printf("\tlog(geff) (p) = %.4f\n",log10(map.leg.eval_00(gsup(),0)(0)));
    printf("\tlog(geff) (e) = %.4f\n",log10(map.leg.eval_00(gsup(),PI/2)(0)));
    printf("\tX0=%.4f   Y0=%.4f   Z0=%.4f\n",X0,Y0,Z0);
    printf("\n");

    printf("Rotation:\n\n");
    double we=map.leg.eval_00(w.row(-1),PI/2)(0);
    double ve=we*re*units.Omega*units.r/1e5;
    printf("\tv_eq = %.3f km/s\n",ve);
    printf("\tOmega (e) = %e rad/s (%.2f%%)\n",we*units.Omega,we/Omegac*100);
    double wp=map.leg.eval_00(w.row(-1),0)(0);
    printf("\tOmega (p) = %e rad/s\n",wp*units.Omega);
    printf("\tOmega (c) = %e rad/s\n",w(0,0)*units.Omega);
    printf("\tPeriod (e) = %.5f days\n",2*PI/we/units.Omega/3600./24.);
    printf("\tPeriod (p) = %.5f days\n",2*PI/wp/units.Omega/3600./24.);
    printf("\tPeriod (c) = %.5f days\n",2*PI/w(0,0)/units.Omega/3600./24.);
    printf("\tLz = %e erg·s\n",Lz());
    printf("\tj  = %e cm2/s (Lz/M)\n",Lz()/M);
    printf("\tT  = %e erg   (Kinetic energy)\n",virial_L()/2);
    printf("\tT/W = %e\n",virial_L()/virial_W()/2);
    printf("\tKelvin-Helmholtz Time = %e yrs\n",virial_3P()/luminosity()/2/365.25/86400.);
    printf("\n");

    if(conv==0) printf("No convective core\n\n");
    else {
        printf("Convective core:\n\n");
        double mcc=Mcore();
        printf("\tMass_core = %.5f Msun (%e g)\n",mcc/M_SUN,mcc);
        double rcc_p=map.leg.eval_00(Rcore(),0)(0);
        printf("\tRadius_core (p) = %.5f Rsun (%e cm)\n",rcc_p/R_SUN,rcc_p);
        double rcc_e=map.leg.eval_00(Rcore(),PI/2)(0);
        printf("\tRadius_core (e) = %.5f Rsun (%e cm)  (flat.=%.3f)\n",rcc_e/R_SUN,rcc_e,1.-rcc_p/rcc_e);
        printf("\tLz_core = %e erg·s\n",Lzcore());
        printf("\tX_core/X_env = %.4f\n",Xc);
        printf("\n");
    }
    printf("Central values:\n\n");
    printf("\tTemperature = %e K\n",Tc);
    printf("\tDensity = %e g/cm3\n",rhoc);
    printf("\tPressure = %e dyn/cm2\n",pc);
    printf("\n");

    printf("Grid parameters:\n\n");
    printf("\t # of domains = %d\n",ndomains);
    printf("\t # of domains in convective core = %d\n",conv);
    printf("\t nr = %d    (",nr);
    for(int n=0;n<ndomains;n++) {
        printf("%d",map.gl.npts[n]);
        if(n<ndomains-1) printf(",");
    }
    printf(")\n");
    printf("\t nth = %d\n",nth);
    printf("\t nex = %d\n",nex);
    printf("\n");

    printf("Additional parameters:\n\n");
    printf("\tOpacity = %s\n",opa.name);
    printf("\tEquation of state = %s\n",eos.name);
    printf("\tNuclear reactions = %s\n",nuc.name);
    printf("\tAtmosphere = %s\n",atm.name);
    printf("\tsurff = %e\n",surff);
    printf("\tcore_convec = %d\n",core_convec);
    printf("\tmin_core_size = %e\n",min_core_size);
    printf("\tenv_convec = %d\n",env_convec);
    printf("\n");

    printf("Tests:\n\n");
    printf("\tVirial test = %e\n",virial());
    printf("\tEnergy test = %e\n",energy_test());

}

void star2d::plot(const matrix_map& error) {

    matrix theta = vector(0, 2*M_PI, 64);
    matrix r = map.leg.eval_00(this->r, theta);
    matrix w = map.leg.eval_00(this->w, theta);

    matrix cost = cos(theta);
    matrix sint = sin(theta);

    matrix x = r*sint;
    matrix y = r*cost;

    plt::clf();

    plt::subplot(231);
    plt::pcolormesh(x, y, w);
    plt::colorbar();
    plt::axis("scaled");

    matrix r_e = map.leg.eval_00(this->r, M_PI/2.0);
    matrix rho_e = map.leg.eval_00(this->rho, M_PI/2.0);
    matrix T_e = map.leg.eval_00(this->T, M_PI/2.0);
    matrix p_e = map.leg.eval_00(this->p, M_PI/2.0);

    plt::subplot(232);
    plt::plot(r_e, rho_e, "$\\rho_{eq}$");
    plt::plot(r_e, T_e, "$T_{eq}$");
    plt::plot(r_e, p_e, "$p_{eq}$");
    plt::legend();


    plt::subplot(233, true);
    std::ostringstream str_stream;

    str_stream.clear();
    str_stream.str("");
    str_stream << Tc;
    plt::text(0.0, .3, std::string("$T_c$:   ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << pc;
    plt::text(0.0, .2, std::string("$p_c$:   ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << pc;
    plt::text(0.0, 0.1, std::string("$\\rho_c$:  ") + str_stream.str());

    str_stream.clear();
    str_stream.str("");
    str_stream << pi_c;
    plt::text(0.0, 0.0, std::string("$\\pi_c$: ") + str_stream.str());

    if (error["Phi"].ncols()*error["Phi"].nrows() > 0) {
        plt::subplot(223);
        plt::semilogy(error["Phi"], "error $\\Phi$");
        plt::semilogy(error["p"], "error $p$");
        plt::legend();
    }

    plt::draw();
    plt::pause();
}

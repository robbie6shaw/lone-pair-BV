CREATE TABLE Ion (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    symbol VARCHAR(2),
    atomic_no INT,
    os INTEGER(1),
    radii REAL,
    softness REAL,
    period INT,
    p_group INT,
    block INT                
);

CREATE TABLE BVParam (
    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
    ion1 INTEGER,
    ion2 INTEGER,
    r0 REAL,
    b REAL,
    ib REAL,
    cn REAL,
    r_cutoff REAL
);
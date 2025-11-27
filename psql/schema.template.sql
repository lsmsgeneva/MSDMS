--
-- PostgreSQL database scheme dump
--

-- Dumped from database version 16.0 (Debian 16.0-1.pgdg120+1)
-- Dumped by pg_dump version 16.0 (Debian 16.0-1.pgdg120+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

SET default_tablespace = '';
SET default_table_access_method = heap;

--
-- Name: metabolites_id_seq; Type: SEQUENCE; Schema: public; Owner: defined in .env
--

CREATE SEQUENCE public.metabolites_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

ALTER SEQUENCE public.metabolites_id_seq OWNER TO {{DB_OWNER}};

--
-- Name: metabolites; Type: TABLE; Schema: public; Owner: defined in .env
--

CREATE TABLE public.metabolites (
                                        id integer DEFAULT nextval('public.metabolites_id_seq'::regclass) NOT NULL PRIMARY KEY,
                                        name character varying,
                                        exact_mass double precision,
                                        chemical_formula character varying,
                                        inchi character varying,
                                        inchikey character varying,
                                        cas_registry_number character varying,
                                        smiles character varying,
                                        pubchem character varying,
                                        chebi character varying,
                                        chemspider character varying
);

CREATE INDEX metabolites_exact_mass_idx ON public.metabolites (exact_mass);

ALTER TABLE public.metabolites OWNER TO {{DB_OWNER}};

--
-- Name: spectra_id_seq; Type: SEQUENCE; Schema: public; Owner: defined in .env
--

CREATE SEQUENCE public.spectra_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

ALTER SEQUENCE public.spectra_id_seq OWNER TO {{DB_OWNER}};

--
-- Name: spectra; Type: TABLE; Schema: public; Owner: defined in .env
--

CREATE TABLE public.spectra (
                                spectrum_id integer DEFAULT nextval('public.spectra_id_seq'::regclass) NOT NULL PRIMARY KEY,
                                metabolite_id integer,
                                splash character varying,
                                instrument_type character varying,
                                collection_date date,
                                ionization_mode character varying,
                                collision_energy_voltage character varying,
                                ms_type character varying,
                                accession character varying,
                                source character varying,
                                fragmentation_method character varying,
                                retention_time double precision,
                                file_path character varying(100),
                                ccs double precision,
                                sv double precision,
                                cov double precision,
                                modifier character varying(5)
);

CREATE INDEX spectra_retention_time_idx ON public.spectra (retention_time);
CREATE INDEX spectra_ccs_idx ON public.spectra (ccs);
CREATE INDEX spectra_sv_idx ON public.spectra (sv);
CREATE INDEX spectra_cov_idx ON public.spectra (cov);
CREATE INDEX spectra_modifier_idx ON public.spectra (modifier);

ALTER TABLE public.spectra OWNER TO {{DB_OWNER}};

--
-- Name: spectrum_peaks; Type: TABLE; Schema: public; Owner: defined in .env
--
CREATE TABLE public.spectrum_peaks (
                                spectrum_id INTEGER REFERENCES public.spectra(spectrum_id) ON DELETE CASCADE,
                                mz_value DOUBLE PRECISION NOT NULL,
                                rel_int DOUBLE PRECISION,
                                PRIMARY KEY (spectrum_id, mz_value)
);

CREATE INDEX idx_mz_value ON public.spectrum_peaks (mz_value);
CREATE INDEX idx_spectrum_id_mz_value ON public.spectrum_peaks (spectrum_id, mz_value);
CREATE INDEX idx_rel_int ON public.spectrum_peaks (rel_int);

ALTER TABLE public.spectrum_peaks OWNER TO {{DB_OWNER}};


--
-- PostgreSQL database dump complete
--


package main

import (
	"compress/gzip"
	"encoding/csv"
	"fmt"
        "io"
	"log"
	"os"
        "strings"
        vk "github.com/genomicsplc/variantkey/go/src"
        "strconv"
        "encoding/gob"
)

func main() {

	args := os.Args[1:]
	gwas_results := args[0]
	filter := args[1]
	false_germline_map := args[2]

        decodeFile, err := os.Open(false_germline_map)
        if err != nil {
            panic(err)
        }
        defer decodeFile.Close()
        decoder := gob.NewDecoder(decodeFile)
        BadKeyStore := make(map[uint64]bool)
        decoder.Decode(&BadKeyStore)
        //filter
	filter_f, err := os.Open(filter)
	if err != nil {
		log.Fatal(err)
	}
	defer filter_f.Close()
	filter_r, err := gzip.NewReader(filter_f)
	if err != nil {
		log.Fatal(err)
	}
	defer filter_r.Close()

	cr := csv.NewReader(filter_r)
        cr.Comma = '\t'

        keys := make(map[uint64]bool)
        count := 0

        for  {

            rec, err := cr.Read()
            if err == io.EOF {
                break
            } else if err != nil {
                log.Fatal(err)
            }

            if count > 0 {

                chrom := rec[0]
                pos, err := strconv.ParseUint(rec[1], 10, 32)
                if err != nil {
                        log.Fatal(err)
                }
                ref := rec[2]
                alt := rec[3]
                key := vk.VariantKey(chrom, uint32(pos), ref, alt)
                keys[key] = true

            }

            count = count + 1
        }

        //gwas
	gwas_f, err := os.Open(gwas_results)
	if err != nil {
		log.Fatal(err)
	}
	defer gwas_f.Close()
	gwas_r, err := gzip.NewReader(gwas_f)
	if err != nil {
		log.Fatal(err)
	}
	defer gwas_r.Close()

	cr = csv.NewReader(gwas_r)
        cr.Comma = '\t'

        count = 0

        for  {

            rec, err := cr.Read()
            if err == io.EOF {
                break
            } else if err != nil {
                log.Fatal(err)
            }

            if count == 0 { // header
                fmt.Println(strings.Join(rec, "\t"))
            }
            if count > 0 {
                // fmt.Println("rec ",  rec)

                chrom := rec[0]
                pos, err := strconv.ParseUint(rec[1], 10, 32)
                if err != nil {
                        log.Fatal(err)
                }
                ref := rec[3]
                alt := rec[4]
                // fmt.Println("len of rec[0] ",  len(rec[0]))
                // fmt.Println("rec[0][:20] ",  rec[0][:5])
                // for _, v := range rec {
                //         fmt.Println("v ", v)
                // }
                key := vk.VariantKey(chrom, uint32(pos), ref, alt)
                if _, ok := keys[key]; ok {
                    // fmt.Println("found record in filter")
                    // fmt.Println("rec ",  rec)
                    continue
                }

                if _, ok := BadKeyStore[key]; ok {
                    // fmt.Println("found record in false germline variant list")
                    // fmt.Println("rec ",  rec)
                    continue
                }
                fmt.Println(strings.Join(rec, "\t"))

            }

            count = count + 1
        }
}

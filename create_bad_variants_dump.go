package main

import (
	"compress/gzip"
	_ "fmt"
        _"io"
	"log"
	"os"
        "strings"
        vk "github.com/genomicsplc/variantkey/go/src"
        _ "strconv"
        "github.com/brentp/vcfgo"
        "path/filepath"
        "math"
        "sync"
        "encoding/gob"
)

type KeyStore struct {
    m map[uint64]bool
    mux sync.Mutex
}

func (store *KeyStore) insert(key uint64) {
    store.mux.Lock()
    // fmt.Printf("now inserting %d\n", key)
    store.m[key] = true
    store.mux.Unlock()
}

func main() {
    dir := "/net/topmed2/working/gt-release/exchange-area/freeze.8/sites/"

    _, err := os.Stat(dir)

    if os.IsNotExist(err) {
        panic(err)
    }

    suffix := "*vcf.gz"

    paths, err := filepath.Glob(filepath.Join(dir, suffix))

    store := KeyStore{m : make(map[uint64]bool)}

    log.Print("before parsing")
    var wg sync.WaitGroup

    for _, path := range paths {
        go parse(path, &store, &wg)
    }

    wg.Wait()

    log.Print("now storing to disk")

    output_file, err := os.Create("likely_somatic_variants_in_freeze8.gob")
    if err != nil {
        panic(err)
    }

    encoder := gob.NewEncoder(output_file)

    if err := encoder.Encode(store.m); err != nil {
        panic(err)
    }

    output_file.Close()

    log.Print("done storing to disk")

}


func parse(path string, store *KeyStore, wg * sync.WaitGroup) {
    wg.Add(1)
    defer wg.Done()
    log.Print("path is ", path)
    f, err := os.Open(path)
    if err != nil {
            log.Fatal(err)
    }
    defer f.Close()
    gzip_file, err := gzip.NewReader(f)
    if err != nil {
            log.Fatal(err)
    }
    defer gzip_file.Close()
    rdr, err := vcfgo.NewReader(gzip_file, false)
    if err != nil {
        panic(err)
    }

    for {
        variant := rdr.Read()

        if variant == nil {
            break
        }

        ABZ_threshold := 5.0 // about 92%
        HWE_SLP_I_threshold := 3.0
        SVM_threshold := -0.25

        alt := strings.Join(variant.Alt(), "")
        key := vk.VariantKey(variant.Chromosome, uint32(variant.Pos), variant.Ref(), alt)
        filter := variant.Filter
        r, err := variant.Info().Get("ABZ")
        ABZ, _ := r.(float64)
        if err != nil {
            panic(err)
        }
        r, err = variant.Info().Get("HWE_SLP_I")
        if err != nil {
            panic(err)
        }
        HWE_SLP_I := r.(float64)
        r, err = variant.Info().Get("SVM")
        if err != nil {
            panic(err)
        }
        SVM := r.(float64)
        //if filter == "PASS" {
        if filter == "PASS" && (math.Abs(ABZ) >= ABZ_threshold || math.Abs(HWE_SLP_I) >= HWE_SLP_I_threshold || SVM < SVM_threshold) {
            // fmt.Printf(
            //     "%d\t%s\t%d\t%s\t%v\t%s\t%.3f\t%.3f\t%.3f\n",
            //     key,
            //     variant.Chromosome,
            //     variant.Pos,
            //     variant.Ref(),
            //     alt,
            //     filter,
            //     ABZ,
            //     HWE_SLP_I,
            //     SVM)
            store.insert(key)
        }
    }
    log.Print("path is now complete ", path)
}

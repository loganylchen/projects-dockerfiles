FROM alpine:3.13.5
VOLUME /circosdata
WORKDIR /circos
# CMD [ "perl", "./your-daemon-or-script.pl" ]
RUN set -x \
    && apk update && apk upgrade \
    && apk add git wget curl libgd perl perl-app-cpanminus tar make gcc build-base gnupg perl-clone perl-extutils-pkgconfig perl-gd perl-devel-checkcompiler perl-module-build-xsutil perl-extutils-config perl-module-build-tiny perl-module-build perl-module-implementation perl-params-validate perl-list-moreutils \
    && cpanm --force Number::Format \
    && cpanm Config::General Font::TTF::Font Math::Bezier Math::Round Math::VecStat Readonly Regexp::Common SVG Set::IntSpan Statistics::Basic Text::Format GD GD::Polyline \
    # && wget http://circos.ca/distribution/circos-tools-0.23.tgz \
    # && tar xfvz circos-tools-0.23.tgz
    && wget http://circos.ca/distribution/circos-0.69-9.tgz \
    && tar xfvz circos-0.69-9.tgz \
    && mv circos-0.69-9/* /circos

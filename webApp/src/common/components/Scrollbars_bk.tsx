import React, { CSSProperties } from 'react';
import Scrollbars from 'react-custom-scrollbars';
import { withTheme } from 'theming';

class CustomScrollbars extends React.Component<any, any> {
  constructor(props, ...rest) {
    super(props, ...rest);
    this.state = {
      scrollTop: 0,
      scrollHeight: 0,
      clientHeight: 0,
    };
  }

  public render() {
    const { scrollTop, scrollHeight, clientHeight } = this.state;
    const shadowTopOpacity = (1 / 20) * Math.min(scrollTop, 20);
    const bottomScrollTop = scrollHeight - clientHeight;
    const shadowBottomOpacity = (1 / 20) * (bottomScrollTop - Math.max(scrollTop, bottomScrollTop - 20));

    const { theme, ...props } = this.props;
    const shadowTopStyle: CSSProperties = {
      position: 'absolute',
      top: 0,
      left: 0,
      right: 0,
      height: 10,
      background: 'linear-gradient(to bottom, rgba(0, 0, 0, 0.2) 0%, rgba(0, 0, 0, 0) 100%)',
      opacity: shadowTopOpacity,
    };
    const shadowBottomStyle: CSSProperties = {
      position: 'absolute',
      bottom: 0,
      left: 0,
      right: 0,
      height: 10,
      background: 'linear-gradient(to top, rgba(0, 0, 0, 0.2) 0%, rgba(0, 0, 0, 0) 100%)',
      opacity: shadowBottomOpacity,
    };
    return (
      <div style={{ position: 'relative', height: '100%' }}>
        <Scrollbars
          renderTrackHorizontal={this.renderTrackHorizontal}
          renderTrackVertical={this.renderTrackVertical}
          renderThumbHorizontal={this.renderThumb}
          renderThumbVertical={this.renderThumb}
          autoHide
          universal
          onUpdate={this.handleUpdate}
          {...props}
        />
        <div style={shadowTopStyle} />
        <div style={shadowBottomStyle} />
      </div>
    );
  }

  private handleUpdate = (values) => {
    const { scrollTop, scrollHeight, clientHeight } = values;
    this.setState({ scrollTop, scrollHeight, clientHeight });
  };

  private renderTrackHorizontal = ({ style, ...props }) => {
    const finalStyle = {
      ...style,
      right: 2,
      bottom: 2,
      left: 2,
      zIndex: 50,
    };
    return <div style={finalStyle} />;
  };

  private renderTrackVertical = ({ style, ...props }) => {
    const finalStyle = {
      ...style,
      right: 2,
      bottom: 2,
      top: 2,
      zIndex: 50,
    };
    return <div style={finalStyle} />;
  };

  private renderThumb = ({ style, ...props }) => {
    const finalStyle = {
      ...style,
      cursor: 'pointer',
      backgroundColor: this.props.theme.palette.text.divider,
    };
    return <div style={finalStyle} />;
  };
}

export default withTheme(CustomScrollbars);
